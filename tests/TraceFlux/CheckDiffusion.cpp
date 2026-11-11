#include "FluxTraceAdvDiv/FluxTraceDiffKernels.hpp"
#include <ADRSolverRedesign/DiffusionScalarIPTraceFlux/DiffusionScalarIPTraceFluxKernels.hpp>
#include <LibUtilities/BasicUtils/Field/MemoryRegion.hpp>
#include <array>
#include <iostream>
#include <stdexcept>

using namespace Nektar;
using namespace PENKNIFE::detail;
using Policy = PENKNIFE::closures::Braginskii<double>;
struct Buffer
{
    LibUtilities::MemoryRegion<double> region;
    double *p;
    explicit Buffer(size_t n) : region(n)
    {
        region.Initialize<NektarSpaces::HostSpace>(0);
        p = region.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
    }
    double &operator[](size_t i) { return p[i]; }
};
void close(double actual, double expected, double magnitude)
{
    if (!std::isfinite(actual) ||
        std::abs(actual - expected) > 1e-10 * std::max(1., magnitude))
        throw std::runtime_error(
            "Anisotropic trace mismatch: " + std::to_string(actual) + " vs " +
            std::to_string(expected));
}

template <typename Exec>
void Check(unsigned dim, const std::array<double, 3> &b,
           const std::array<double, 3> &normal)
{
    constexpr size_t N = 12;
    // Sparse global field indices; electron energy is deliberately first.
    constexpr unsigned comps = 7, ni = 1, mi = 3, ui = 4, ue = 0;
    unsigned nIds[] = {ni}, vIds[] = {mi}, eIds[] = {ui};
    double mass[] = {4}, charge[] = {2};
    TraceDiffusionData<Policy, double> data{{1e18, 100, constants::qeomp, 1},
                                            1,
                                            ue,
                                            nIds,
                                            vIds,
                                            eIds,
                                            mass,
                                            charge,
                                            .4,
                                            .7,
                                            1.1,
                                            2.3,
                                            1.7};
    Buffer mean(comps * N), left(comps * N), right(comps * N),
        calibrationB(3 * N), calibrationQ(comps * dim * N),
        calibrationFlux(comps * dim * N), B(3 * N), normals(dim * N),
        penalty(N), qf(comps * dim * N), qb(comps * dim * N), nff(comps * N),
        nfb(comps * N), weights(comps * N), out(comps * N), wsState(comps * N),
        wsGradient(comps * dim * N), wsFlux(comps * dim * N);
    TraceDiffusionWorkspace<double> workspace{wsState.p, wsGradient.p,
                                              wsFlux.p};
    for (size_t i = 0; i < N - 1; ++i)
    {
        mean[ni * N + i] = 1.5 + .03 * i;
        // Signed momentum must not change the density-based diffusivity.
        mean[mi * N + i] = (i % 2 == 0 ? 1. : -1.) * (2. + .02 * i);
        mean[ui * N + i] = 4.5 + .05 * i;
        mean[ue * N + i] = 9. + .04 * i;
        for (unsigned f = 0; f < comps; ++f)
        {
            const double jump = f == ni   ? .2 + .01 * i
                                : f == mi ? -.3
                                : f == ui ? .8
                                : f == ue ? -.4
                                          : 0;
            left[f * N + i] = mean[f * N + i] - .5 * jump;
            right[f * N + i] = mean[f * N + i] + .5 * jump;
            weights[f * N + i] = f == ui ? 0 : f == mi ? .5 : 1;
            nff[f * N + i] = .3 * (f + 1) + .01 * i;
            nfb[f * N + i] = -.2 * (f + 1) - .02 * i;
            for (unsigned d = 0; d < dim; ++d)
            {
                qf[(f * dim + d) * N + i] = .04 * (f + 1) * (d + 1);
                qb[(f * dim + d) * N + i] = -.02 * (f + 2) * (d + 1);
            }
        }
        calibrationB[i] = 1;
        calibrationQ[(ui * dim) * N + i] = 1;
        calibrationQ[(ue * dim) * N + i] = 1;
        for (unsigned d = 0; d < 3; ++d)
            B[d * N + i] = b[d];
        for (unsigned d = 0; d < dim; ++d)
            normals[d * N + i] = normal[d];
        penalty[i] = 1.2 + .01 * i;
    }
    // Obtain each mean-state parallel coefficient from the existing closure
    // with unit energy gradients and no density gradients. This avoids
    // duplicating its collision-frequency model in the test.
    FluxVolDiffKernel<false, Exec, Policy>(
        data.closureParameters, N, dim, comps, 1, ue, nIds, vIds, eIds, mass,
        charge, 0., 0., 0., N, N, N, N, mean.p, calibrationQ.p,
        calibrationFlux.p, calibrationB.p, 1., 0);

    for (unsigned form = 0; form < 3; ++form)
    {
        for (size_t i = 0; i < comps * N; ++i)
            out[i] = .9;
        if (form == 0)
            DiffuseTraceFluxKernel<Exec, Policy, true>(
                data, N, dim, comps, normals.p, penalty.p,
                static_cast<const double *>(nullptr), B.p, left.p, right.p,
                qf.p, qb.p, out.p, workspace);
        else if (form == 1)
            AvgNormalFluxTraceKernel<Exec, Policy>(
                data, N, dim, comps, normals.p, penalty.p, B.p, left.p, right.p,
                nff.p, nfb.p, out.p, workspace);
        else
            DiffuseTraceFluxKernel<Exec, Policy, false>(
                data, N, dim, comps, normals.p, penalty.p, weights.p, B.p,
                left.p, right.p, qf.p, qb.p, out.p, workspace);

        for (size_t i = 0; i < N - 1; ++i)
        {
            const double ki = calibrationFlux[(ui * dim) * N + i];
            const double ke = calibrationFlux[(ue * dim) * N + i];
            if (!(ki > 0 && ke > 0 && std::abs(ki - ke) > 1e-12))
                throw std::runtime_error(
                    "Expected distinct ion/electron heat coefficients");
            for (unsigned f = 0; f < comps; ++f)
            {
                double reference = 0, magnitude = 1;
                // Apply the independent ADR tensor kernel to every nonzero
                // block D_fj. Summing its scalar answers gives the coupled
                // plasma numerical flux, including cross-field penalties.
                for (unsigned j = 0; j < comps; ++j)
                {
                    double parallel = 0, perpendicular = 0;
                    if (f == ni && j == ni)
                        perpendicular = data.anomD;
                    if (f == mi && j == mi)
                        perpendicular = data.anomNu * mean[ni * N + i];
                    if (f == mi && j == ni)
                        perpendicular = data.anomD;
                    if (f == ui && j == ui)
                    {
                        parallel = ki;
                        perpendicular = data.anomChi * mean[ui * N + i];
                    }
                    if (f == ui && j == ni)
                    {
                        parallel = -ki * mean[ui * N + i] / mean[ni * N + i];
                        perpendicular = data.anomD;
                    }
                    if (f == ue && j == ue)
                        parallel = ke;
                    if (f == ue && j == ni)
                        parallel = -ke * mean[ue * N + i] / mean[ni * N + i];
                    if (parallel == 0 && perpendicular == 0)
                        continue;
                    double tensor[6] = {}, g0[3] = {}, g1[3] = {}, value = 0;
                    for (unsigned d = 0; d < dim; ++d)
                    {
                        if (form != 1)
                        {
                            g0[d] = qf[(j * dim + d) * N + i];
                            g1[d] = qb[(j * dim + d) * N + i];
                        }
                        for (unsigned e = 0; e < dim; ++e)
                            tensor[LibUtilities::GetDiffCoeffMap(dim,
                                                                 d * dim + e)] =
                                data.scale *
                                (parallel * b[d] * b[e] +
                                 perpendicular *
                                     ((d == e ? 1. : 0.) - b[d] * b[e]));
                    }
                    const double l = left[j * N + i], r = right[j * N + i];
                    Nektar::detail::DiffuseScalarTraceFluxKernel<
                        NektarSpaces::Serial, true>(
                        1, dim, 1, 1, data.penaltyCoeff, tensor, normal.data(),
                        penalty.p + i, &l, &r, g0, g1, &value);
                    reference += value;
                    magnitude += std::abs(value);
                }
                if (form == 1)
                    reference += .5 * (nff[f * N + i] - nfb[f * N + i]);
                if (form == 2)
                    reference *= weights[f * N + i];
                close(out[f * N + i], -.9 + reference, magnitude);
            }
        }
        // Zero padding must neither evaluate log(0) nor create any diffusion.
        for (unsigned f = 0; f < comps; ++f)
            close(out[f * N + N - 1], -.9, 1);
    }
}

template <typename Exec> void Run()
{
    Check<Exec>(1, {1, 0, 0}, {1, 0, 0});
    Check<Exec>(1, {.6, 0, .8}, {1, 0, 0});
    Check<Exec>(2, {1, 0, 0}, {1, 0, 0});
    Check<Exec>(2, {1, 0, 0}, {0, 1, 0});
    Check<Exec>(2, {.6, .8, 0}, {.8, .6, 0});
    Check<Exec>(2, {.6, 0, .8}, {.6, .8, 0});
    Check<Exec>(3, {0, 0, 1}, {0, 0, 1});
    Check<Exec>(3, {0, 0, 1}, {1, 0, 0});
    Check<Exec>(3, {.6, 0, .8}, {0, .6, .8});
}

int main()
{
    Run<NektarSpaces::Serial>();
#if defined(__AVX2__)
    Run<NektarSpaces::AVX>();
    std::cout << "Serial and AVX2 diffusion checks passed.\n";
#else
    std::cout << "Serial diffusion checks passed.\n";
#endif
    std::cout << "Compared with ADR tensor blocks in 1D/2D/3D: parallel, "
                 "perpendicular, oblique, field-specific strengths, density "
                 "coupling, both forms, boundaries and padding.\n";
}
