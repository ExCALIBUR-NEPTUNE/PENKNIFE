#pragma once

#include "LibUtilities/LoopExecution/LoopExecution.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename TData>
NEK_FORCE_INLINE static void FluxVolAdvKernel(
    const unsigned int npts, const unsigned int n_species,
    const unsigned int *n_ids, const unsigned int *v_ids,
    const unsigned int *e_ids, const TData *masses, const unsigned int velComps,
    const unsigned int BStride, const unsigned int velStride,
    const unsigned int inStride, const unsigned int outStride,
    const TData *velbase, const TData *inbase, TData *outbase,
    const TData *Bbase, const TData scale, const unsigned int streamID)
{
    using vec_t =
        typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                              TData>::type;
    constexpr unsigned int vec_width =
        (std::is_same_v<ExecSpace, NektarSpaces::AVX>)
            ? tinysimd::simd<TData>::width
            : 1;

    const size_t groupsize    = npts / vec_width;
    const size_t velVecStride = velStride / vec_width;
    const size_t inVecStride  = inStride / vec_width;
    const size_t outVecStride = outStride / vec_width;
    const size_t BVecStride   = BStride / vec_width;

    const auto velvec = reinterpret_cast<const vec_t *>(velbase);
    const auto invec  = reinterpret_cast<const vec_t *>(inbase);
    auto outvec       = reinterpret_cast<vec_t *>(outbase);
    const auto Bvec   = reinterpret_cast<const vec_t *>(Bbase);

    Nektar::LoopExecutionSetStreamID(streamID);

    // Parallelize over point groups; each i is independent.
    Nektar::parallel_for<ExecSpace>(
        0u, groupsize, NEKTAR_LAMBDA(const size_t i) {
            vec_t b_vec[3] = {vec_t(TData(0)), vec_t(TData(0)),
                              vec_t(TData(0))};
            for (unsigned int d = 0; d < ndim; ++d)
            {
                b_vec[d] = Bvec[3 * BVecStride + d];
            }
            for (unsigned int spec = 0; spec < n_species; ++spec)
            {
                const unsigned int n_idx = n_ids[spec];
                const unsigned int v_idx = v_ids[spec];
                const unsigned int e_idx = e_ids[spec];
                const TData mass         = masses[spec];

                const size_t n_off = n_idx * inVecStride + i;
                const vec_t n_vec  = invec[n_off];

                const size_t v_off  = v_idx * inVecStride + i;
                const vec_t mom_vec = invec[v_off];
                const vec_t v_vec   = SafeDivide(mom_vec, n_vec * mass);

                const size_t e_off = e_idx * inVecStride + i;
                const vec_t e_vec  = invec[e_off];
                const vec_t p_vec  = (2. / 3) * e_vec;

                for (unsigned int d = 0; d < velComps; ++d)
                {
                    const size_t vel_off =
                        spec * velComps * velVecStride + d * velVecStride + i;
                    const vec_t vp_vec = velvec[vel_off];
                    const size_t out_off =
                        (n_idx * velComps + d) * outVecStride + i;

                    outvec[out_off] +=
                        scale * n_vec * (v_vec * b_vec[d] + vp_vec);

                    out_off = (v_idx * velComps + d) * outVecStride + i;
                    outvec[out_off] +=
                        scale * (mom_vec * (v_vec * b_vec[d] + vp_vec) +
                                 p_vec * b_vec[d]);

                    out_off = (e_idx * velComps + d) * outVecStride + i;
                    outvec[out_off] +=
                        scale * e_val * (v_vec * b_vec[d] + vp_vec);
                }
            }
        });

    Nektar::LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail