#pragma once

#include <LibUtilities/LoopExecution/LoopExecution.hpp>

namespace PENKNIFE::detail
{
using namespace Nektar;

template <template <typename> typename FluxKernel, typename ExecSpace,
          typename TData>
NEK_FORCE_INLINE static void
FluxKernelLauncher(const size_t npts, const unsigned ndim,
                   const unsigned nvarComps, const unsigned nSpecies,
                   const unsigned electronEnergyId, const unsigned *nIds,
                   const unsigned *vIds, const unsigned *eIds,
                   const TData *masses, const TData *charges,
                   const TData *normal, const TData *b, const TData *currentFwd,
                   const TData *currentBwd, const TData *fwd, const TData *bwd,
                   TData *flux, const TData scale, const unsigned streamID = 0)
{
    constexpr unsigned width = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                   ? tinysimd::simd<TData>::width
                                   : 1;
    ASSERTL0(ndim >= 1 && ndim <= 3 && electronEnergyId < nvarComps,
             "Invalid advection trace dimension or electron-energy index.");
    ASSERTL1(npts % width == 0,
             "The padded trace block must be a whole number of vectors.");
    LoopExecutionSetStreamID(streamID);
    parallel_for<ExecSpace>(
        0u, npts / width, NEKTAR_LAMBDA(const size_t i) {
            const size_t offset = i * width;
            FluxKernel<ExecSpace>()(
                npts, ndim, nvarComps, nSpecies, electronEnergyId, nIds, vIds,
                eIds, masses, charges, normal + offset, b + offset,
                currentFwd + offset, currentBwd + offset, fwd + offset,
                bwd + offset, flux + offset, scale);
        });
    LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail
