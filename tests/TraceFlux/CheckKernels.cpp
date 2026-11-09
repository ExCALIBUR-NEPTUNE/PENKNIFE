#include "FluxTraceAdvDiv/FluxTraceAdvDiffOp.hpp"
#include "FluxTraceAdvDiv/FluxTraceAdvKernels.hpp"
#include "FluxTraceAdvDiv/FluxTraceDiffKernels.hpp"
#include "FluxVolAdvDiff/FluxVolAdvDiffKernels.hpp"
#include "RiemannSolvers/LaxFriedrichsSolverKernels.hpp"
#include <LibUtilities/BasicUtils/Field/MemoryRegion.hpp>
#include <iostream>
#include <stdexcept>

using namespace Nektar;
using namespace PENKNIFE::detail;
using Policy = PENKNIFE::closures::Braginskii<double>;
using Exec = NektarSpaces::Serial;
void close(double a, double b)
{
    if (!std::isfinite(a) ||
        std::abs(a - b) > 1e-10 * std::max(1., std::abs(b)))
        throw std::runtime_error("Mismatch: " + std::to_string(a) + " vs " +
                                 std::to_string(b));
}
struct Buffer
{
    LibUtilities::MemoryRegion<double> region;
    double *p;
    Buffer(size_t n) : region(n)
    {
        region.Initialize<NektarSpaces::HostSpace>(0);
        p = region.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
    }
    double &operator[](size_t i) { return p[i]; }
};
int main()
{
    constexpr size_t N = 12;
    constexpr unsigned comps = 10, dim = 2, ns = 3, ee = 0;
    constexpr unsigned bad = std::numeric_limits<unsigned>::max();
    unsigned ni[] = {1, 3, 5}, vi[] = {2, bad, 6}, ei[] = {4, bad, 7};
    double masses[] = {1, 2, 4}, charges[] = {1, 0, 2};
    Buffer left(comps * N), right(comps * N), b(3 * N), normal(dim * N), j0(N),
        j1(N), flux(comps * N), rev(comps * N), ne(N), me(N),
        q(comps * dim * N), fv(comps * dim * N), total(comps * dim * N);
    for (size_t i = 0; i < N - 1; ++i)
    {
        left[0 * N + i] = 36;
        left[1 * N + i] = 2;
        left[2 * N + i] = 3;
        left[3 * N + i] = 1e5;
        left[4 * N + i] = 3;
        left[5 * N + i] = 3;
        left[6 * N + i] = -4;
        left[7 * N + i] = 9;
        b[i] = .6;
        b[2 * N + i] = .8;
        normal[i] = .8;
        normal[N + i] = .6;
        j0[i] = j1[i] = .2 + .1 * i;
        ne[i] = 8;
        me[i] = constants::m_e_m_p * (3. - 2. - j0[i]);
    }
    std::copy(left.p, left.p + comps * N, right.p);
    auto launch = [&](double *f, double *r, double *out, double *jf, double *jr,
                      double scale = 1.)
    {
        FluxKernelLauncher<LaxFriedrichsSolverKernel, Exec>(
            N, dim, comps, ns, ee, ni, vi, ei, masses, charges, normal.p, b.p,
            jf, jr, f, r, out, scale);
    };
    launch(left.p, right.p, flux.p, j0.p, j1.p, 2.);
    FluxVolAdvDiffKernel<Exec, Policy>(
        {1e18, 100, constants::qeomp, 1}, N, dim, comps, ns, ee, ni, vi, ei,
        masses, charges, constants::m_e_m_p, 0., 0., 0., N, N, N, N, N, left.p,
        q.p, fv.p, total.p, b.p, ne.p, me.p, 2., 0);
    for (unsigned f = 0; f < comps; ++f)
        for (size_t i = 0; i < N; ++i)
        {
            double expected = 0;
            for (unsigned d = 0; d < dim; ++d)
                expected -= normal[d * N + i] * total[(f * dim + d) * N + i];
            close(flux[f * N + i], expected);
        }
    // Unequal states: conservative under side exchange + normal reversal.
    for (size_t i = 0; i < N - 1; ++i)
    {
        for (unsigned f = 0; f < comps; ++f)
            right[f * N + i] *= 1.3;
        right[2 * N + i] = -2;
        right[6 * N + i] = 5;
        j1[i] = -.4;
    }
    launch(left.p, right.p, flux.p, j0.p, j1.p);
    for (size_t i = 0; i < dim * N; ++i)
        normal[i] = -normal[i];
    launch(right.p, left.p, rev.p, j1.p, j0.p);
    for (size_t i = 0; i < comps * N; ++i)
        close(flux[i], -rev[i]);
    // Tangential normal: both flux and numerical viscosity must vanish.
    for (size_t i = 0; i < N; ++i)
    {
        normal[i] = 0;
        normal[N + i] = 1;
    }
    launch(left.p, right.p, flux.p, j0.p, j1.p);
    for (size_t i = 0; i < comps * N; ++i)
        close(flux[i], 0);
    // Prescribed current changes electron advection only for equal states.
    for (size_t i = 0; i < N; ++i)
    {
        normal[i] = 1;
        normal[N + i] = 0;
    }
    launch(left.p, left.p, flux.p, j0.p, j0.p);
    for (size_t i = 0; i < N - 1; ++i)
        j1[i] = j0[i] + 2;
    launch(left.p, left.p, rev.p, j1.p, j1.p);
    for (unsigned f = 0; f < comps; ++f)
        for (size_t i = 0; i < N - 1; ++i)
            close(rev[f * N + i] - flux[f * N + i],
                  f == ee ? -left[ee * N + i] * .6 * 2 / ne[i] : 0);
    // All-zero padded lanes stay finite and zero.
    for (unsigned f = 0; f < comps; ++f)
        close(rev[f * N + N - 1], 0);
    // Runtime registration uses the same public key as Create().
    if (!MultiRegions::GetOperatorFactory<double>().ModuleExists(
            "FluxTraceAdvDiffLaxFriedrichsBraginskiiSerial"))
        throw std::runtime_error("Trace factory key missing");
    // Boundary ghost helper must preserve pre-reflected states.
    Buffer g0(8), g1(8), gd(8);
    unsigned flags[] = {0, 1};
    for (unsigned i = 0; i < 8; ++i)
    {
        g0[i] = 1;
        g1[i] = 3;
        gd[i] = 3;
    }
    ReflectDiffusionGhostStateKernel<Exec>(2, 1, 2, 4, flags, g0.p, g1.p, gd.p);
    for (unsigned f = 0; f < 4; ++f)
    {
        close(gd[2 * f], 5);
        close(gd[2 * f + 1], 3);
    }
    std::cout << "Trace/advection checks passed: volume consistency, current, "
                 "reversal, tangential limit, padding, factory and boundary "
                 "ghosts.\n";
}
