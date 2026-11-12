#pragma once
#include <LibUtilities/BasicUtils/SafeDivide.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <limits>

namespace PENKNIFE::detail
{
using namespace Nektar;

template <typename ExecSpace, typename TData>
NEK_FORCE_INLINE static void FluxVolAdvKernel(
    const size_t npts, const unsigned nDim, const unsigned nvarComps,
    const unsigned n_species, const unsigned electronEnergyId,
    const unsigned *n_ids, const unsigned *v_ids, const unsigned *e_ids,
    const TData *masses, const TData *charges, const TData electronMass,
    const size_t inStride, const size_t outStride, const size_t BStride,
    const size_t electronStride, const TData *inbase, TData *outbase,
    const TData *Bbase, const TData *nebase, const TData *Mebase,
    const TData scale, const unsigned streamID)
{
    using vec_t =
        typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                              TData>::type;
    constexpr unsigned width     = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                       ? tinysimd::simd<TData>::width
                                       : 1;
    constexpr unsigned invalidId = std::numeric_limits<unsigned>::max();
    ASSERTL0(nDim >= 1 && nDim <= 3 && electronEnergyId < nvarComps &&
                 electronMass > 0,
             "Invalid dimensions, electron energy index or electron mass.");
    ASSERTL1(npts <= inStride && npts <= derivStride && npts <= outStride &&
                 npts <= BStride && npts <= electronStride &&
                 npts % width == 0 && inStride % width == 0 &&
                 outStride % width == 0 && BStride % width == 0 &&
                 electronStride % width == 0,
             "Combined flux block strides must match the SIMD width.");
    const auto invec = reinterpret_cast<const vec_t *>(inbase);
    const auto bvec  = reinterpret_cast<const vec_t *>(Bbase);
    const auto nevec = reinterpret_cast<const vec_t *>(nebase);
    const auto Mevec = reinterpret_cast<const vec_t *>(Mebase);
    auto outvec      = reinterpret_cast<vec_t *>(outbase);

    const size_t inputStride  = inStride / width;
    const size_t outputStride = outStride / width;
    const size_t bStride      = BStride / width;
    LoopExecutionSetStreamID(streamID);
    parallel_for<ExecSpace>(
        0u, npts / width, NEKTAR_LAMBDA(const size_t i) {
            for (unsigned s = 0; s < n_species; ++s)
            {
                if (n_ids[s] == invalidId)
                    continue;
                const vec_t n        = invec[n_ids[s] * inputStride + i];
                const vec_t momentum = v_ids[s] == invalidId
                                           ? vec_t(0)
                                           : invec[v_ids[s] * inputStride + i];
                const vec_t energy   = e_ids[s] == invalidId
                                           ? vec_t(0)
                                           : invec[e_ids[s] * inputStride + i];
                const vec_t velocity =
                    SafeDivide(momentum, vec_t(masses[s]) * n);
                const vec_t pressure = (TData(2) / 3) * energy;
                for (unsigned d = 0; d < nDim; ++d)
                {
                    const vec_t b      = bvec[d * bStride + i];
                    const vec_t advVel = velocity * b;
                    outvec[(n_ids[s] * nDim + d) * outputStride + i] -=
                        scale * n * advVel;
                    if (v_ids[s] != invalidId)
                        outvec[(v_ids[s] * nDim + d) * outputStride + i] -=
                            scale * (momentum * advVel + pressure * b);
                    if (e_ids[s] != invalidId)
                        outvec[(e_ids[s] * nDim + d) * outputStride + i] -=
                            scale * energy * advVel;
                }
            }
            const vec_t ue =
                SafeDivide(Mevec[i], vec_t(electronMass) * nevec[i]);
            const vec_t Ue = invec[electronEnergyId * inputStride + i];
            for (unsigned d = 0; d < nDim; ++d)
                outvec[(electronEnergyId * nDim + d) * outputStride + i] -=
                    scale * Ue * ue * bvec[d * bStride + i];
        });
    LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail
