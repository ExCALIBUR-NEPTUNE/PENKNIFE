#include "BndCond/BndCondBohm/BndCondBohmKernels.hpp"
#include "FluxTraceAdvDiv/FluxTraceDiffKernels.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
using namespace Nektar;
using namespace PENKNIFE;
using namespace PENKNIFE::detail;
struct Buffer
{
    LibUtilities::MemoryRegion<double> region;
    double *p;
    Buffer(size_t n) : region(n)
    {
        region.Initialize<NektarSpaces::HostSpace>(0);
        p = region.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
    }
    double &operator[](size_t i)
    {
        return p[i];
    }
};
void close(double a, double b)
{
    if (!std::isfinite(a) ||
        std::abs(a - b) > 2e-10 * std::max({1., std::abs(a), std::abs(b)}))
        throw std::runtime_error("Mismatch " + std::to_string(a) + " vs " +
                                 std::to_string(b));
}
void check(unsigned dim, std::array<double, 3> b, std::array<double, 3> normal,
           double u, bool secondIon, double current = 0,
           double secondCharge = 2)
{
    constexpr size_t N   = 4;
    constexpr unsigned C = 8, ee = 0;
    unsigned ns[] = {1, 5}, ms[] = {2, 6}, es[] = {3, 7};
    double masses[] = {2, 4}, charges[] = {1, secondCharge};
    const unsigned count = secondIon ? 2 : 1;
    BohmSpecies<double> species{count, ee, ns, ms, es, masses, charges, 1.7};
    Buffer state(C * N), normals(dim * N), auxiliary(4 * N);
    double bn = 0, ne = 0;
    for (unsigned d = 0; d < dim; ++d)
        bn += b[d] * normal[d];
    for (unsigned s = 0; s < count; ++s)
        ne += charges[s] * (1.3 + s);
    for (size_t i = 0; i < N; ++i)
    {
        for (unsigned d = 0; d < dim; ++d)
            normals[d * N + i] = normal[d];
        for (unsigned d = 0; d < 3; ++d)
            auxiliary[d * N + i] = b[d];
        auxiliary[3 * N + i] = current;
        if (i == N - 1)
            continue;
        for (unsigned s = 0; s < count; ++s)
        {
            state[ns[s] * N + i] = 1.3 + s;
            state[ms[s] * N + i] = masses[s] * (1.3 + s) * u;
            state[es[s] * N + i] = 1.5 * (1.3 + s) * (2 + s);
        }
        state[ee * N + i] = 1.5 * ne * 3.;
        state[4 * N + i]  = 42.;
    }
    ApplyBohmBoundary<NektarSpaces::Serial>(species, BohmParameters<double>{},
                                            N, dim, normals.p, auxiliary.p,
                                            state.p);
    for (size_t i = 0; i < N - 1; ++i)
    {
        double electrons = -current * bn;
        for (unsigned s = 0; s < count; ++s)
        {
            const double n = 1.3 + s, Ti = 2 + s;
            const double cs = std::sqrt((charges[s] * 3. + 3 * Ti) / masses[s]);
            const double sign = bn > 0 ? 1 : -1;
            const double ub =
                charges[s] != 0 && bn != 0 ? sign * std::max(sign * u, cs) : u;
            const double particles = n * ub * bn;
            close(state[ns[s] * N + i], 1.7 * particles);
            close(state[ms[s] * N + i],
                  1.7 * (masses[s] * n * ub * ub + n * Ti) * bn);
            close(state[es[s] * N + i],
                  1.7 * (charges[s] != 0 ? 2.5 : 1.5) * n * Ti * ub * bn);
            electrons += charges[s] * particles;
        }
        close(state[ee * N + i], 1.7 * 4.5 * 3. * electrons);
        close(state[4 * N + i], 42.); // unrelated field is left alone
    }
    for (unsigned f = 0; f < C; ++f)
        close(state[f * N + N - 1], 0);
}
// A prescribed complete boundary flux must not call the transport closure.
struct RejectClosure
{
    struct Parameters
    {
    };
    template <typename... Args> static void Evaluate(Args...)
    {
        throw std::runtime_error("Prescribed sheath flux evaluated a closure");
    }
};
void checkNoTransport()
{
    constexpr size_t N   = 4;
    constexpr unsigned C = 4;
    Buffer weights(C * N), flux(C * N);
    std::fill(flux.p, flux.p + C * N, 7.);
    TraceDiffusionData<RejectClosure, double> data{};
    TraceDiffusionWorkspace<double> workspace{};
    const double *unused = nullptr;
    DiffuseTraceFluxKernel<NektarSpaces::Serial, RejectClosure, false>(
        data, N, 1, C, unused, unused, weights.p, unused, unused, unused,
        unused, unused, flux.p, workspace);
    for (size_t i = 0; i < C * N; ++i)
        close(flux[i], -7.);
}
int main()
{
    try
    {
        checkNoTransport();
        for (double u : {-12., -1., 0., 1., 12.})
        {
            check(1, {1, 0, 0}, {1, 0, 0}, u, false);
            check(2, {-.6, .8, 0}, {1, 0, 0}, u, true);
            check(3, {.6, 0, .8}, {.6, .8, 0}, u, true, .2);
            check(2, {.3, .4, std::sqrt(.75)}, {.6, .8, 0}, u, true);
            check(2, {0, .6, .8}, {1, 0, 0}, u, true);
            check(2, {.6, .8, 0}, {1, 0, 0}, u, true, .2, 0);
        }
        std::cout << "Bohm speed, free density, total sheath energy loss, "
                     "current, neutrals and oblique B checks passed.\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
