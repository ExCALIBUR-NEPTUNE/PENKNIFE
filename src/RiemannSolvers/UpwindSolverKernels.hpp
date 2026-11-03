#pragma once

// The dimension and shape kernels. NOTE: They are NOT duplicate
// templated version based on the array size like the
// operators. HOWEVER, they are forced to be INLINED. The inlining is
// critical so that when used in the templated version of the operator
// that loop unrolling occurs.

namespace PENKNIFE::detail
{

template <typename ExecSpace> struct UpwindSolverKernel
{

    template <typename TScalar>
    NEK_DEVICE_INLINE void operator()(
        const size_t blksize, const unsigned int n_species,
        const unsigned int *n_ids, const unsigned int *v_ids,
        const unsigned int *e_ids, const TScalar *masses,
        const TScalar *charges, const unsigned velComps, const TScalar *velbase,
        const TScalar *normbase, const TScalar *fwdBbase,
        const TScalar *bwdBbase, const TScalar *fwdnebase,
        const TScalar *bwdnebase, const TScalar *fwdbase,
        const TScalar *bwdbase, TScalar *fluxbase)
    {
        using std::abs;

        // Explicit vectorisation for AVX backend, vec_t =
        // tinysimd::simd<TScalar> for AVX, vec_t = TScalar otherwise.
        using vec_t =
            typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                                  TScalar>::type;
        const unsigned vec_width =
            (std::is_same_v<ExecSpace, NektarSpaces::AVX>)
                ? tinysimd::simd<TScalar>::width
                : 1;

        // Assumes blksize is a whole number of vectors; checked on the host
        // by FluxKernelLauncher, since a kernel may not assert.

        // Parallelize over points; each i is independent
        const size_t groupsize = blksize / vec_width;
        const auto velvec      = reinterpret_cast<const vec_t *>(velbase);
        const auto normvec     = reinterpret_cast<const vec_t *>(normbase);
        const auto fwdvec      = reinterpret_cast<const vec_t *>(fwdbase);
        const auto bwdvec      = reinterpret_cast<const vec_t *>(bwdbase);
        auto fluxvec           = reinterpret_cast<vec_t *>(fluxbase);

        for (unsigned int spec = 0; spec < n_species; ++spec)
        {
            const unsigned int n_idx = n_ids[spec];
            const unsigned int v_idx = v_ids[spec];
            const unsigned int e_idx = e_ids[spec];
            // Build nv = v·n on the fly
            vec_t nv = 0.0;
            for (unsigned d = 0; d < velComps; ++d)
            {
                const size_t vel_off =
                    spec * velComps * groupsize + d * groupsize;
                const vec_t v_i = velvec[vel_off];
                const vec_t n_i = normvec[d * groupsize];
                nv += v_i * n_i;
            }

            const vec_t nv_abs = abs(nv);

            // Branchless split: nv_pos=max(nv,0), nv_neg=min(nv,0)
            const vec_t nv_pos = 0.5 * (nv + nv_abs);
            const vec_t nv_neg = 0.5 * (nv - nv_abs);

            // Blend without branches:
            // flux = nv_pos * Fwd + nv_neg * Bwd

            fluxvec[n_idx * groupsize] = nv_pos * fwdvec[n_idx * groupsize] +
                                         nv_neg * bwdvec[n_idx * groupsize];
            fluxvec[v_idx * groupsize] = nv_pos * fwdvec[v_idx * groupsize] +
                                         nv_neg * bwdvec[v_idx * groupsize];
            fluxvec[e_idx * groupsize] = nv_pos * fwdvec[e_idx * groupsize] +
                                         nv_neg * bwdvec[e_idx * groupsize];
        }
    }
};

} // namespace PENKNIFE::detail
