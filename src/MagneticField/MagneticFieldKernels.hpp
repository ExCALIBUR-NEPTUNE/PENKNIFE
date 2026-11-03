#pragma once

#include <LibUtilities/BasicUtils/SafeDivide.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <cmath>
#include <type_traits>

namespace PENKNIFE::detail
{
using namespace Nektar;
template <typename ExecSpace, typename TData>
NEK_FORCE_INLINE void MagneticFieldKernel(
    const size_t npts, const size_t BStride, const size_t unitStride,
    const TData *B, TData *unit, TData *magSquared, const unsigned int streamID)
{
    using vec_t = typename Nektar::data_type_if<
        std::is_same_v<ExecSpace, NektarSpaces::AVX>, TData>::type;
    constexpr size_t width = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                 ? tinysimd::simd<TData>::width
                                 : 1;
    ASSERTL1(npts % width == 0 && BStride % width == 0 &&
                 unitStride % width == 0,
             "Magnetic field blocks must be padded to the vector width.");
    const auto in   = reinterpret_cast<const vec_t *>(B);
    auto out        = reinterpret_cast<vec_t *>(unit);
    auto mag        = reinterpret_cast<vec_t *>(magSquared);
    const size_t bs = BStride / width;
    const size_t us = unitStride / width;

    Nektar::LoopExecutionSetStreamID(streamID);
    Nektar::parallel_for<ExecSpace>(
        0u, npts / width, NEKTAR_LAMBDA(const size_t i) {
            using std::sqrt;
            const vec_t bx      = in[i];
            const vec_t by      = in[bs + i];
            const vec_t bz      = in[2 * bs + i];
            const vec_t squared = bx * bx + by * by + bz * bz;
            const vec_t norm    = sqrt(squared);
            // Overwrite on every update; also keep zero padding finite.
            mag[i]          = squared;
            out[i]          = Nektar::SafeDivide(bx, norm);
            out[us + i]     = Nektar::SafeDivide(by, norm);
            out[2 * us + i] = Nektar::SafeDivide(bz, norm);
        });
    Nektar::LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail
