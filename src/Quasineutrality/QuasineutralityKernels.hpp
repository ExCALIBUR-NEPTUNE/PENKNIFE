#pragma once

#include "LibUtilities/LoopExecution/LoopExecution.hpp"

namespace PENKNIFE::detail
{
using namespace Nektar;
template <typename ExecSpace, typename TData>
NEK_FORCE_INLINE static void QuasineutralKernel(
    const size_t npts, const unsigned int n_species, const unsigned int *n_ids,
    const unsigned int *momentumIds, const TData *charges,
    const TData *momentumWeights, const unsigned int invalidId,
    const size_t inStride, const size_t outStride, const TData electronMass,
    const TData *inbase, const TData *currentBase, TData *densityBase,
    TData *momentumBase, const unsigned int streamID)
{
    using vec_t =
        typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                              TData>::type;
    constexpr unsigned int vec_width =
        (std::is_same_v<ExecSpace, NektarSpaces::AVX>)
            ? tinysimd::simd<TData>::width
            : 1;

    ASSERTL1(npts <= inStride && npts <= outStride && npts % vec_width == 0 &&
                 inStride % vec_width == 0,
             "Quasineutrality point counts and strides must match the vector "
             "width.");
    const size_t groupsize   = npts / vec_width;
    const size_t inVecStride = inStride / vec_width;

    const auto invec      = reinterpret_cast<const vec_t *>(inbase);
    const auto currentVec = reinterpret_cast<const vec_t *>(currentBase);
    auto densityVec       = reinterpret_cast<vec_t *>(densityBase);
    auto momentumVec      = reinterpret_cast<vec_t *>(momentumBase);

    Nektar::LoopExecutionSetStreamID(streamID);

    // Parallelize over point groups; each i is independent.
    Nektar::parallel_for<ExecSpace>(
        0u, groupsize, NEKTAR_LAMBDA(const size_t i) {
            vec_t density  = vec_t(0.0);
            vec_t momentum = vec_t(0.0);
            for (unsigned int spec = 0; spec < n_species; ++spec)
            {
                const TData charge       = charges[spec];
                const unsigned int n_idx = n_ids[spec];

                const size_t n_off = n_idx * inVecStride + i;
                const vec_t n_vec  = invec[n_off];
                density += charge * n_vec;
                const unsigned int momentumId = momentumIds[spec];
                if (momentumId != invalidId)
                {
                    momentum += momentumWeights[spec] *
                                invec[momentumId * inVecStride + i];
                }
            }
            densityVec[i]  = density;
            momentumVec[i] = momentum - electronMass * currentVec[i];
        });

    Nektar::LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail