#pragma once

#include "FluxVolDiff/FluxVolDiffKernels.hpp"

namespace PENKNIFE::detail
{
// Same state ordering, normalization and transport parameters as the volume
// flux. All metadata pointers must be accessible in ExecSpace.
template <typename ClosureType, typename TData> struct TraceDiffusionData
{
    typename ClosureType::Parameters closureParameters;
    unsigned nSpecies, electronEnergyId;
    const unsigned *nIds, *vIds, *eIds;
    const TData *masses, *charges;
    TData anomD, anomNu, anomChi;
    TData penaltyCoeff = 1;
    TData scale = 1;
};

// Caller-owned, aligned, non-aliasing arrays sized for the padded block:
// state: nvarComps*npts; gradient and flux: nvarComps*ndim*npts each.
// Each iteration owns its point/SIMD group's entries in all three arrays.
// The point flux evaluator reads only those entries, so no inter-thread
// synchronization is needed inside the fused launch. Reuse this storage
// between blocks on the same operator stream.
template <typename TData> struct TraceDiffusionWorkspace
{
    TData *state, *gradient, *flux;
};

/**
 * Assemble -F_adv + F_diff.n using the volume diffusion law F(U,q,b).
 *
 * Ubar = (Uf + Ub)/2, jump = Ub - Uf, tau = penaltyCoeff*penaltyFactor.
 * AverageNormalFluxes:
 *   Fhat_diff = (nFf - nFb)/2 + n.F(Ubar, tau*n*jump, b).
 * Otherwise:
 *   Fhat_diff = n.F(Ubar, (qf+qb)/2 + tau*n*jump, b).
 *
 * nFb uses the backward element's outward normal, hence the minus sign.
 * The magnetic unit vector is continuous across the trace. Use its forward
 * trace directly; do not average or renormalize its reduced-dimensional
 * components (a 2D mesh may still have an out-of-plane magnetic component).
 *
 * F is linear in the gradients at fixed state, including the temperature
 * gradients ((2/3)*grad(U_s)-T_s*grad(n_s))/n_s and the quasineutral electron
 * density gradient. Applying it to the rank-one penalty gradient therefore
 * retains the full coupled anisotropic penalty rather than a scalar penalty
 * applied independently to each conserved field.
 *
 * This is the system analogue of the redesigned ADR solver's
 * tau*(n^T D n)*[u] penalty. For F_f = sum_j D_fj grad(U_j), it is
 * tau*sum_j (n^T D_fj n)*[U_j]. The blocks D_fj depend on the field:
 * density uses anomalous D, momentum uses anomalous nu*nbar, and ion
 * internal energy combines anomalous chi*Ubar with parallel conduction.
 * Electron internal energy uses its own parallel conductivity. Density
 * gradients also couple into momentum and heat conduction.
 * Writing Ppar=(b.n)^2 and Pperp=|n|^2-Ppar, perpendicular blocks acquire
 * Pperp and parallel blocks acquire Ppar. For example the density penalty
 * is tau*anomD*Pperp*[n_s]; an ion's parallel heat penalty is
 * tau*kpar_s*Ppar*((2/3)*[U_s]-Tbar_s*[n_s])/nbar_s.
 * The electron term uses [n_e]=sum_s Z_s*[n_s]. These coefficients are
 * evaluated at the mean state, as selected for both trace formulations.
 * No additional common scalar penalty or side-maximum coefficient is
 * applied: that would change the selected law or count diffusion twice.
 *
 * Boundary callers provide the diffusion ghost state and gradients, just as
 * in the redesigned CFE solver. This gives the same half-average and Ub-Uf
 * jump on physical boundaries as on interior and partition traces.
 * Boundary weights are component-major (nvarComps*npts); they affect only
 * diffusion. Interior kernels do not read them. Boundary traces always use
 * the gradient-based formulation, even when interior traces average fluxes.
 */
template <bool AverageNormalFluxes, bool IsInterior, typename ExecSpace,
          typename ClosureType, typename TData>
NEK_FORCE_INLINE static void DiffusionTraceKernel(
    const TraceDiffusionData<ClosureType, TData> data, const size_t npts,
    const unsigned ndim, const unsigned nvarComps, const TData *normbase,
    const TData *penaltyFactorBase, const TData *boundaryFluxWeightBase,
    const TData *Bbase, const TData *fwdbase, const TData *bwdbase,
    const TData *secondfwdbase, const TData *secondbwdbase, TData *fluxbase,
    const TraceDiffusionWorkspace<TData> workspace, const unsigned streamID)
{
    static_assert(!AverageNormalFluxes || IsInterior,
                  "Boundary diffusion requires gradient traces.");
    using vec_t =
        typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                              TData>::type;
    constexpr unsigned width = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                   ? tinysimd::simd<TData>::width
                                   : 1;
    ASSERTL0(ndim >= 1 && ndim <= 3 && data.electronEnergyId < nvarComps,
             "Invalid trace dimension or electron-energy index.");
    ASSERTL0(IsInterior || boundaryFluxWeightBase,
             "Boundary diffusion requires a weight for every component.");
    ASSERTL1(npts % width == 0,
             "The padded trace block must be a whole number of vectors.");
    const size_t stride = npts / width;
    const auto fwd = reinterpret_cast<const vec_t *>(fwdbase);
    const auto bwd = reinterpret_cast<const vec_t *>(bwdbase);
    const auto normal = reinterpret_cast<const vec_t *>(normbase);
    const auto penalty = reinterpret_cast<const vec_t *>(penaltyFactorBase);
    const auto secondFwd = reinterpret_cast<const vec_t *>(secondfwdbase);
    const auto secondBwd = reinterpret_cast<const vec_t *>(secondbwdbase);
    auto state = reinterpret_cast<vec_t *>(workspace.state);
    auto gradient = reinterpret_cast<vec_t *>(workspace.gradient);

    const auto physicalFlux = reinterpret_cast<const vec_t *>(workspace.flux);
    const auto weights =
        reinterpret_cast<const vec_t *>(boundaryFluxWeightBase);
    auto flux = reinterpret_cast<vec_t *>(fluxbase);

    LoopExecutionSetStreamID(streamID);
    parallel_for<ExecSpace>(
        0u, stride, NEKTAR_LAMBDA(const size_t i) {
            const vec_t tau = data.penaltyCoeff * penalty[i];
            for (unsigned f = 0; f < nvarComps; ++f)
            {
                const size_t offset = f * stride + i;
                // Energies are already internal-energy densities. Do not
                // reconstruct CFE total energy or assume ndim+2 variables.
                state[offset] = TData(0.5) * (fwd[offset] + bwd[offset]);
                const vec_t jump = bwd[offset] - fwd[offset];
                for (unsigned d = 0; d < ndim; ++d)
                {
                    const size_t qoffset = (f * ndim + d) * stride + i;
                    vec_t q = tau * normal[d * stride + i] * jump;
                    if constexpr (!AverageNormalFluxes)
                        q += TData(0.5) *
                             (secondFwd[qoffset] + secondBwd[qoffset]);
                    gradient[qoffset] = q;
                }
            }

            // State, gradients, closure evaluation and normal contraction
            // all belong to this point/SIMD group. Call the point evaluator
            // directly to avoid launching a separate volume-flux kernel.
            FluxVolDiffPoint<false, width, ClosureType, vec_t>(
                data.closureParameters, i, ndim, nvarComps, data.nSpecies,
                data.electronEnergyId, data.nIds, data.vIds, data.eIds,
                data.masses, data.charges, data.anomD, data.anomNu,
                data.anomChi, npts, npts, npts, npts, workspace.state,
                workspace.gradient, workspace.flux, Bbase, data.scale);

            for (unsigned f = 0; f < nvarComps; ++f)
            {
                const size_t offset = f * stride + i;
                vec_t diffusive = vec_t(0);
                for (unsigned d = 0; d < ndim; ++d)
                    diffusive += normal[d * stride + i] *
                                 physicalFlux[(f * ndim + d) * stride + i];
                if constexpr (AverageNormalFluxes)
                    // Input normal fluxes already carry the volume scale.
                    diffusive +=
                        TData(0.5) * (secondFwd[offset] - secondBwd[offset]);
                if constexpr (!IsInterior)
                    diffusive *= weights[offset];
                flux[offset] = -flux[offset] + diffusive;
            }
        });
    LoopExecutionSetStreamID(0);
}

template <typename ExecSpace, typename ClosureType, typename TData>
NEK_FORCE_INLINE static void AvgNormalFluxTraceKernel(
    const TraceDiffusionData<ClosureType, TData> data, const size_t npts,
    const unsigned ndim, const unsigned nvarComps, const TData *normbase,
    const TData *penaltyFactorBase, const TData *Bbase, const TData *fwdbase,
    const TData *bwdbase, const TData *nFfwdbase, const TData *nFbwdbase,
    TData *fluxbase, const TraceDiffusionWorkspace<TData> workspace,
    const unsigned streamID = 0)
{
    DiffusionTraceKernel<true, true, ExecSpace>(
        data, npts, ndim, nvarComps, normbase, penaltyFactorBase,
        static_cast<const TData *>(nullptr), Bbase, fwdbase, bwdbase, nFfwdbase,
        nFbwdbase, fluxbase, workspace, streamID);
}

template <typename ExecSpace, typename ClosureType, bool IsInterior,
          typename TData>
NEK_FORCE_INLINE static void DiffuseTraceFluxKernel(
    const TraceDiffusionData<ClosureType, TData> data, const size_t npts,
    const unsigned ndim, const unsigned nvarComps, const TData *normbase,
    const TData *penaltyFactorBase, const TData *boundaryFluxWeightBase,
    const TData *Bbase, const TData *fwdbase, const TData *bwdbase,
    const TData *derivfwdbase, const TData *derivbwdbase, TData *fluxbase,
    const TraceDiffusionWorkspace<TData> workspace, const unsigned streamID = 0)
{
    DiffusionTraceKernel<false, IsInterior, ExecSpace>(
        data, npts, ndim, nvarComps, normbase, penaltyFactorBase,
        boundaryFluxWeightBase, Bbase, fwdbase, bwdbase, derivfwdbase,
        derivbwdbase, fluxbase, workspace, streamID);
}
// Preserve already-reflected ghosts; reflect imposed states only.
template <typename ExecSpace, typename TData>
NEK_FORCE_INLINE static void ReflectDiffusionGhostStateKernel(
    const size_t numTrace, const size_t npTot, const size_t npTBlock,
    const unsigned numflux, const unsigned *ghostIsReflected,
    const TData *gloT0, const TData *gloT1, TData *gloT1Diff,
    const unsigned int streamID = 0)
{
    Nektar::LoopExecutionSetStreamID(streamID);

    const size_t npTrace = numTrace * npTot;

    Nektar::parallel_for<ExecSpace>(
        0u, npTrace * numflux, NEKTAR_LAMBDA(const size_t idx) {
            const size_t f = idx / npTrace;
            const size_t r = idx - f * npTrace;
            const size_t t = r / npTot;

            if (ghostIsReflected[t])
            {
                return;
            }

            const size_t o = f * npTBlock + r;

            gloT1Diff[o] = TData(2.0) * gloT1[o] - gloT0[o];
        });

    Nektar::LoopExecutionSetStreamID(0);
}

} // namespace PENKNIFE::detail
