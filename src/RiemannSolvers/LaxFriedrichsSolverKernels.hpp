#pragma once

#include <LibUtilities/BasicUtils/SafeDivide.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace PENKNIFE::detail
{
using namespace Nektar;

// Local Lax-Friedrichs flux for the same internal-energy state used by the
// combined volume flux. Input pointers address one SIMD group; stride is
// still the full padded block length. J_parallel is in number-current units,
// exactly as in QuasineutralKernel: ne*ue = sum_s Z_s*M_s/m_s - J_parallel.
template <typename ExecSpace> struct LaxFriedrichsSolverKernel
{
    template <typename TData>
    NEK_HOSTDEVICE_INLINE void
    operator()(const size_t npts, const unsigned ndim, const unsigned nvarComps,
               const unsigned nSpecies, const unsigned electronEnergyId,
               const unsigned *nIds, const unsigned *vIds, const unsigned *eIds,
               const TData *masses, const TData *charges,
               const TData *normalBase, const TData *Bbase,
               const TData *currentFwdBase, const TData *currentBwdBase,
               const TData *fwdBase, const TData *bwdBase, TData *fluxBase,
               const TData scale) const
    {
        using std::abs;
        using std::max;
        using std::sqrt;
        using vec_t =
            typename data_type_if<std::is_same_v<ExecSpace, NektarSpaces::AVX>,
                                  TData>::type;
        constexpr unsigned width = std::is_same_v<ExecSpace, NektarSpaces::AVX>
                                       ? tinysimd::simd<TData>::width
                                       : 1;
        constexpr unsigned invalid = std::numeric_limits<unsigned>::max();
        const size_t stride = npts / width;
        const auto normal = reinterpret_cast<const vec_t *>(normalBase);
        const auto b = reinterpret_cast<const vec_t *>(Bbase);
        const auto fwd = reinterpret_cast<const vec_t *>(fwdBase);
        const auto bwd = reinterpret_cast<const vec_t *>(bwdBase);
        auto flux = reinterpret_cast<vec_t *>(fluxBase);
        const vec_t half = vec_t(TData(0.5));
        const vec_t twoThirds = vec_t(TData(2) / 3);
        vec_t bn = vec_t(0), neL = vec_t(0), neR = vec_t(0);
        vec_t particleL = vec_t(0), particleR = vec_t(0);
        for (unsigned d = 0; d < ndim; ++d)
            bn += b[d * stride] * normal[d * stride];
        for (unsigned s = 0; s < nSpecies; ++s)
        {
            if (charges[s] == 0 || nIds[s] == invalid)
                continue;
            neL += charges[s] * fwd[nIds[s] * stride];
            neR += charges[s] * bwd[nIds[s] * stride];
            if (vIds[s] != invalid)
            {
                particleL += (charges[s] / masses[s]) * fwd[vIds[s] * stride];
                particleR += (charges[s] / masses[s]) * bwd[vIds[s] * stride];
            }
        }
        const vec_t UeL = fwd[electronEnergyId * stride];
        const vec_t UeR = bwd[electronEnergyId * stride];
        const vec_t TeL = SafeDivide(twoThirds * UeL, neL);
        const vec_t TeR = SafeDivide(twoThirds * UeR, neR);
        const vec_t ueL =
            bn * SafeDivide(particleL - *reinterpret_cast<const vec_t *>(
                                            currentFwdBase),
                            neL);
        const vec_t ueR =
            bn * SafeDivide(particleR - *reinterpret_cast<const vec_t *>(
                                            currentBwdBase),
                            neR);
        vec_t electronSpeed = max(abs(ueL), abs(ueR));
        for (unsigned f = 0; f < nvarComps; ++f)
            flux[f * stride] = vec_t(0);
        for (unsigned s = 0; s < nSpecies; ++s)
        {
            if (nIds[s] == invalid)
                continue;
            const unsigned ni = nIds[s], vi = vIds[s], ei = eIds[s];
            const vec_t nL = fwd[ni * stride], nR = bwd[ni * stride];
            const vec_t ML = vi == invalid ? vec_t(0) : fwd[vi * stride];
            const vec_t MR = vi == invalid ? vec_t(0) : bwd[vi * stride];
            const vec_t UL = ei == invalid ? vec_t(0) : fwd[ei * stride];
            const vec_t UR = ei == invalid ? vec_t(0) : bwd[ei * stride];
            const vec_t uL = bn * SafeDivide(ML, vec_t(masses[s]) * nL);
            const vec_t uR = bn * SafeDivide(MR, vec_t(masses[s]) * nR);
            const vec_t pL = twoThirds * UL, pR = twoThirds * UR;
            const vec_t TiL = SafeDivide(pL, nL), TiR = SafeDivide(pR, nR);
            // Retain the existing ion-acoustic model, with its normal
            // projection. Max side spectral radii define the local LF speed.
            const vec_t cL = sqrt(
                max(vec_t(0), (charges[s] * TeL + TData(3) * TiL) / masses[s]));
            const vec_t cR = sqrt(
                max(vec_t(0), (charges[s] * TeR + TData(3) * TiR) / masses[s]));
            const vec_t a = max(abs(uL) + abs(bn) * cL, abs(uR) + abs(bn) * cR);
            electronSpeed = max(electronSpeed, a);
            flux[ni * stride] =
                scale * half * (nL * uL + nR * uR - a * (nR - nL));
            if (vi != invalid)
                flux[vi * stride] =
                    scale * half *
                    (ML * uL + pL * bn + MR * uR + pR * bn - a * (MR - ML));
            if (ei != invalid)
                flux[ei * stride] =
                    scale * half * (UL * uL + UR * uR - a * (UR - UL));
        }
        flux[electronEnergyId * stride] =
            scale * half *
            (UeL * ueL + UeR * ueR - electronSpeed * (UeR - UeL));
    }
};
} // namespace PENKNIFE::detail
