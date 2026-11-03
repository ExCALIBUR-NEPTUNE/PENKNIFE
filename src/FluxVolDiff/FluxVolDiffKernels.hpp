#pragma once

#include "Closures/BraginskiiHeatFlux.hpp"
#include <LibUtilities/BasicUtils/SafeDivide.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <limits>

namespace PENKNIFE::detail
{
using namespace Nektar;

template <bool APPEND, unsigned vec_width, typename ClosureType, typename vec_t,
          typename TData>
NEK_HOSTDEVICE_INLINE static void FluxVolDiffPoint(
    const typename ClosureType::Parameters &closureParameters, const size_t i,
    const unsigned nDim, const unsigned nvarComps, const unsigned n_species,
    const unsigned electronEnergyId, const unsigned *n_ids,
    const unsigned *v_ids, const unsigned *e_ids, const TData *masses,
    const TData *charges, const TData anom_D, const TData anom_nu,
    const TData anom_chi, const size_t inStride, const size_t derivStride,
    const size_t outStride, const size_t BStride, const TData *inbase,
    const TData *qbase, TData *outbase, const TData *Bbase, const TData scale)
{
    constexpr unsigned invalidId = std::numeric_limits<unsigned>::max();
    // Species arrays use the global vector ordering. Charged species require
    // density and energy; uncharged species do not enter the closure sums.
    const auto invec            = reinterpret_cast<const vec_t *>(inbase);
    const auto qvec             = reinterpret_cast<const vec_t *>(qbase);
    auto outvec                 = reinterpret_cast<vec_t *>(outbase);
    const auto Bvec             = reinterpret_cast<const vec_t *>(Bbase);
    const size_t inVecStride    = inStride / vec_width;
    const size_t derivVecStride = derivStride / vec_width;
    const size_t outVecStride   = outStride / vec_width;
    const size_t BVecStride     = BStride / vec_width;

    if constexpr (!APPEND)
    {
        for (unsigned f = 0; f < nvarComps * nDim; ++f)
            outvec[f * outVecStride + i] = vec_t(0);
    }
    // Closure: collision frequencies and parallel ion/electron heat
    // fluxes. Evaluate each SIMD lane independently for the Coulomb-log
    // branches. All workspace is point-local, including frequencies and
    // gradients.
    for (unsigned lane = 0; lane < vec_width; ++lane)
        ClosureType::Evaluate(closureParameters, i * vec_width + lane, nDim,
                              n_species, electronEnergyId, n_ids, e_ids, masses,
                              charges, inStride, derivStride, outStride,
                              BStride, inbase, qbase, Bbase, outbase, scale);

    // Anomalous diffusions follow the closure and add to its heat
    // fluxes. Keep the existing anomalous expressions; repair their
    // component and direction indices here without changing their
    // transport model.
    for (unsigned spec = 0; spec < n_species; ++spec)
    {
        const unsigned n_idx = n_ids[spec];
        const unsigned v_idx = v_ids[spec];
        const unsigned e_idx = e_ids[spec];
        if (n_idx == invalidId)
            continue;
        for (unsigned outDir = 0; outDir < nDim; ++outDir)
        {
            for (unsigned derivDir = 0; derivDir < nDim; ++derivDir)
            {
                const vec_t bb = vec_t(outDir == derivDir ? 1 : 0) -
                                 Bvec[outDir * BVecStride + i] *
                                     Bvec[derivDir * BVecStride + i];
                const size_t n_derivOff =
                    (n_idx * nDim + derivDir) * derivVecStride + i;
                outvec[(n_idx * nDim + outDir) * outVecStride + i] +=
                    bb * scale * anom_D * qvec[n_derivOff];
                if (v_idx != invalidId)
                {
                    const size_t v_derivOff =
                        (v_idx * nDim + derivDir) * derivVecStride + i;
                    outvec[(v_idx * nDim + outDir) * outVecStride + i] +=
                        bb * scale *
                        (anom_nu * invec[n_idx * inVecStride + i] *
                             qvec[v_derivOff] +
                         anom_D * qvec[n_derivOff]);
                }
                if (e_idx != invalidId)
                {
                    const size_t e_derivOff =
                        (e_idx * nDim + derivDir) * derivVecStride + i;
                    outvec[(e_idx * nDim + outDir) * outVecStride + i] +=
                        bb * scale *
                        (anom_chi * invec[e_idx * inVecStride + i] *
                             qvec[e_derivOff] +
                         anom_D * qvec[n_derivOff]);
                }
            }
        }
    }
}

template <bool APPEND, typename ExecSpace, typename ClosureType, typename TData>
NEK_FORCE_INLINE static void FluxVolDiffKernel(
    const typename ClosureType::Parameters closureParameters, const size_t npts,
    const unsigned int nDim, const unsigned int nvarComps,
    const unsigned int n_species, const unsigned int electronEnergyId,
    const unsigned int *n_ids, const unsigned int *v_ids,
    const unsigned int *e_ids, const TData *masses, const TData *charges,
    const TData anom_D, const TData anom_nu, const TData anom_chi,
    const size_t inStride, const size_t derivStride, const size_t outStride,
    const size_t BStride, const TData *inbase, const TData *qbase,
    TData *outbase, const TData *Bbase, const TData scale,
    const unsigned int streamID)
{
    using vec_t =
        typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                              TData>::type;
    constexpr unsigned vec_width = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                       ? tinysimd::simd<TData>::width
                                       : 1;
    ASSERTL0(nDim >= 1 && nDim <= 3 && electronEnergyId < nvarComps,
             "Invalid dimensions or electron-energy index for diffusion flux.");
    ASSERTL1(npts <= inStride && npts <= derivStride && npts <= outStride &&
                 npts <= BStride && npts % vec_width == 0 &&
                 inStride % vec_width == 0 && derivStride % vec_width == 0 &&
                 outStride % vec_width == 0 && BStride % vec_width == 0,
             "Diffusion block strides must match the SIMD width.");
    LoopExecutionSetStreamID(streamID);
    parallel_for<ExecSpace>(
        0u, npts / vec_width, NEKTAR_LAMBDA(const size_t i) {
            FluxVolDiffPoint<APPEND, vec_width, ClosureType, vec_t>(
                closureParameters, i, nDim, nvarComps, n_species,
                electronEnergyId, n_ids, v_ids, e_ids, masses, charges, anom_D,
                anom_nu, anom_chi, inStride, derivStride, outStride, BStride,
                inbase, qbase, outbase, Bbase, scale);
        });
    LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail
