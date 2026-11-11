#pragma once

#include "BndCondBohmOp.hpp"
#include <LibUtilities/BasicUtils/SafeDivide.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <cmath>
#include <limits>

namespace PENKNIFE::detail
{
template <typename T> struct BohmParameters
{
    T acousticGamma   = 3;
    T ionTransmission = T(2.5), electronTransmission = T(4.5);
};

// Transform the seeded state in place, like the CFS boundary block kernels.
// Density is extrapolated; temperatures come from internal energy. Neutrals
// retain extrapolated advection. Only mapped species/electron slots are used.
template <typename ExecSpace, typename T>
void ApplyBohmBoundary(BohmSpecies<T> species, BohmParameters<T> parameters,
                       size_t stride, unsigned dim, const T *normal,
                       const T *magneticCurrent, T *boundary,
                       unsigned streamID = 0)
{
    LoopExecutionSetStreamID(streamID);
    parallel_for<ExecSpace>(
        0u, stride, NEKTAR_LAMBDA(size_t i) {
            using std::sqrt;
            constexpr unsigned invalid = std::numeric_limits<unsigned>::max();
            T bn = 0, ne = 0;
            for (unsigned d = 0; d < dim; ++d)
                bn += normal[d * stride + i] * magneticCurrent[d * stride + i];
            for (unsigned s = 0; s < species.count; ++s)
                if (species.charge[s] != 0)
                    ne += species.charge[s] *
                          boundary[species.density[s] * stride + i];
            const unsigned ee = species.electronEnergy;
            const T Te  = SafeDivide(T(2) / 3 * boundary[ee * stride + i], ne);
            T electrons = -magneticCurrent[3 * stride + i] * bn;
            for (unsigned s = 0; s < species.count; ++s)
            {
                const unsigned ni = species.density[s],
                               mi = species.momentum[s], ei = species.energy[s];
                if (ni == invalid)
                    continue;
                const T n = boundary[ni * stride + i];
                const T M = mi == invalid ? 0 : boundary[mi * stride + i];
                const T U = ei == invalid ? 0 : boundary[ei * stride + i];
                const T p = T(2) / 3 * U;
                T u       = SafeDivide(M, species.mass[s] * n);
                if (species.charge[s] != 0 && n != 0 && bn != 0)
                {
                    const T cs   = sqrt((species.charge[s] * Te +
                                       parameters.acousticGamma * p / n) /
                                        species.mass[s]);
                    const T sign = bn > 0 ? 1 : -1;
                    u            = sign * (sign * u > cs ? sign * u : cs);
                }
                const T particles         = n * u * bn;
                boundary[ni * stride + i] = species.scale * particles;
                if (mi != invalid)
                    boundary[mi * stride + i] =
                        species.scale * (species.mass[s] * n * u * u + p) * bn;
                if (ei != invalid)
                    boundary[ei * stride + i] =
                        species.scale *
                        (species.charge[s] != 0 ? parameters.ionTransmission * p
                                                : U) *
                        u * bn;
                electrons += species.charge[s] * particles;
            }
            // Replace total thermal-energy flux, including existing advection.
            boundary[ee * stride + i] = species.scale *
                                        parameters.electronTransmission * Te *
                                        electrons;
        });
    LoopExecutionSetStreamID(0);
}
} // namespace PENKNIFE::detail
