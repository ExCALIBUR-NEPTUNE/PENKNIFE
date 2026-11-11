#pragma once

#include "BraginskiiHeatFlux.hpp"
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <limits>

namespace PENKNIFE::closures
{
// The scalar models use the same volume/trace evaluator. B is the physical
// magnetic field for these two policies, so its magnitude is retained.
template <typename T> struct SingleDiffusive
{
    static inline const std::string name = "SingleDiffusive";
    static constexpr bool RequiresEnergy = false;
    struct Parameters
    {
        T Nnorm, Tnorm, omega_c, mesh_length;
        T parallel, perpendicular;
        bool fixedParallel, fixedPerpendicular;
    };
    static void Load(const Nektar::LibUtilities::SessionReaderSharedPtr &s,
                     Parameters &p)
    {
        T kc, lambda, temperature;
        s->LoadParameter("k_c", kc);
        s->LoadParameter("lambda", lambda);
        s->LoadParameter("T_bg", temperature);
        p.parallel = 6 * kc * sqrt(2 * pow(M_PI, 3)) / lambda *
                     constants::epsilon_0 * constants::epsilon_0 *
                     constants::c * pow(p.Tnorm * temperature, 2.5) /
                     sqrt(constants::m_e) * 1e12 / p.Nnorm;
        p.perpendicular =
            lambda / (6 * sqrt(p.Tnorm * temperature * pow(M_PI, 3))) *
            sqrt(constants::m_p) / constants::c *
            pow(constants::e / constants::epsilon_0_si, 2) * p.Nnorm;
        const T scale = p.omega_c * p.mesh_length * p.mesh_length;
        p.parallel /= scale;
        p.perpendicular /= scale;
        p.fixedParallel      = s->DefinesParameter("k_par");
        p.fixedPerpendicular = s->DefinesParameter("k_perp");
        if (p.fixedParallel)
            s->LoadParameter("k_par", p.parallel);
        if (p.fixedPerpendicular)
            s->LoadParameter("k_perp", p.perpendicular);
    }
    NEK_HOSTDEVICE_INLINE static void Evaluate(
        const Parameters &p, size_t point, unsigned dim, unsigned ns, unsigned,
        const unsigned *nIds, const unsigned *, const T *masses,
        const T *charges, size_t is, size_t qs, size_t os, size_t bs,
        const T *in, const T *q, const T *B, T *out, T scale)
    {
        T b[3], B2 = 0;
        for (unsigned d = 0; d < 3; ++d)
        {
            b[d] = B[d * bs + point];
            B2 += b[d] * b[d];
        }
        if (B2 == 0)
            return;
        for (auto &x : b)
            x /= sqrt(B2);
        for (unsigned s = 0; s < ns; ++s)
        {
            if (charges[s] == 0 ||
                nIds[s] == std::numeric_limits<unsigned>::max())
                continue;
            const T n = in[nIds[s] * is + point];
            if (n <= 0)
                continue;
            const T par  = p.fixedParallel
                               ? p.parallel
                               : p.parallel / (charges[s] * charges[s] * n);
            const T perp = p.fixedPerpendicular
                               ? p.perpendicular
                               : p.perpendicular * charges[s] * charges[s] *
                                     sqrt(masses[s]) * n / B2;
            for (unsigned d = 0; d < dim; ++d)
                for (unsigned j = 0; j < dim; ++j)
                    out[(nIds[s] * dim + d) * os + point] +=
                        scale *
                        ((par - perp) * b[d] * b[j] + (d == j ? perp : 0)) *
                        q[(nIds[s] * dim + j) * qs + point];
        }
    }
};

template <typename T> struct DoubleDiffusive
{
    static inline const std::string name = "DoubleDiffusive";
    static constexpr bool RequiresEnergy = true;
    struct Parameters
    {
        T Nnorm, Tnorm, omega_c, mesh_length;
        T kpar, kperp, kipar, kiperp, kicross, kepar, keperp, kecross, kci;
    };
    static void Load(const Nektar::LibUtilities::SessionReaderSharedPtr &s,
                     Parameters &p)
    {
        T lambda, kce, Bnorm;
        s->LoadParameter("lambda", lambda);
        s->LoadParameter("k_ci", p.kci, 3.9);
        s->LoadParameter("k_ce", kce, 3.16);
        s->LoadParameter("Bnorm", Bnorm, 1.0);
        const T scaling = p.omega_c * p.mesh_length * p.mesh_length;
        const T tau     = 6 * sqrt(2 * pow(M_PI, 3)) / lambda *
                      constants::epsilon_0 * constants::epsilon_0 * 1e12 *
                      pow(p.Tnorm, 1.5) / (constants::c * p.Nnorm);
        const T gyro = 1 / (constants::c * constants::c);
        p.kpar = kce * 1.5 * p.Tnorm * p.Nnorm * tau / sqrt(constants::m_e) /
                 scaling;
        p.kperp = 1.5 * p.Tnorm * p.Nnorm * gyro * gyro / tau / scaling;
        p.kipar = p.kci * 1.5 * p.Tnorm * p.Nnorm * tau / sqrt(constants::m_p) /
                  scaling;
        p.kiperp = 3 * p.Tnorm * p.Nnorm * gyro * gyro * sqrt(constants::m_p) /
                   tau / scaling;
        p.kicross = 3.75 * p.Tnorm * p.Nnorm * gyro / Bnorm / scaling;
        p.kepar   = kce * 1.5 * p.Tnorm * p.Nnorm * tau / sqrt(constants::m_e) /
                  scaling;
        p.keperp = (sqrt(2.) + 3.25) * 1.5 * p.Tnorm * p.Nnorm * gyro * gyro *
                   sqrt(constants::m_e) / tau / scaling;
        p.kecross = p.kicross;
    }
    NEK_HOSTDEVICE_INLINE static void Evaluate(
        const Parameters &p, size_t point, unsigned dim, unsigned ns,
        unsigned ee, const unsigned *nIds, const unsigned *eIds,
        const T *masses, const T *charges, size_t is, size_t qs, size_t os,
        size_t bs, const T *in, const T *q, const T *B, T *out, T scale)
    {
        constexpr unsigned invalid = std::numeric_limits<unsigned>::max();
        T b[3], B2 = 0, ne = 0, gradNe[3] = {0, 0, 0};
        for (unsigned d = 0; d < 3; ++d)
        {
            b[d] = B[d * bs + point];
            B2 += b[d] * b[d];
        }
        if (B2 == 0)
            return;
        for (auto &x : b)
            x /= sqrt(B2);
        for (unsigned s = 0; s < ns; ++s)
            if (charges[s] != 0 && nIds[s] != invalid)
            {
                ne += charges[s] * in[nIds[s] * is + point];
                for (unsigned d = 0; d < dim; ++d)
                    gradNe[d] +=
                        charges[s] * q[(nIds[s] * dim + d) * qs + point];
            }
        if (ne <= 0)
            return;
        const T Te = (T(2) / 3) * in[ee * is + point] / ne;
        if (Te <= 0)
            return;
        for (unsigned s = 0; s <= ns; ++s)
        {
            const bool electron = s == ns;
            if (!electron &&
                (charges[s] == 0 || nIds[s] == invalid || eIds[s] == invalid))
                continue;
            const unsigned energy = electron ? ee : eIds[s];
            const T n             = electron ? ne : in[nIds[s] * is + point];
            const T temperature   = (T(2) / 3) * in[energy * is + point] / n;
            if (n <= 0 || temperature <= 0)
                continue;
            T par, perp, cross;
            if (electron)
            {
                par   = p.kepar * pow(Te, 2.5);
                perp  = p.keperp * ne * ne / (B2 * sqrt(Te));
                cross = p.kecross * ne * Te / sqrt(B2);
            }
            else
            {
                const T Z = charges[s], A = masses[s];
                const T dp = p.kci * p.kpar * pow(Te, 2.5) / (Z * Z * n);
                const T dt = p.kperp * Z * Z * sqrt(A) * n / (sqrt(Te) * B2);
                for (unsigned d = 0; d < dim; ++d)
                    for (unsigned j = 0; j < dim; ++j)
                        out[(nIds[s] * dim + d) * os + point] +=
                            scale *
                            ((dp - dt) * b[d] * b[j] + (d == j ? dt : 0)) *
                            q[(nIds[s] * dim + j) * qs + point];
                T sum = 0;
                for (unsigned r = 0; r < ns; ++r)
                    if (charges[r] != 0 && nIds[r] != invalid)
                        sum += charges[r] * charges[r] *
                               sqrt(masses[r] / (A + masses[r])) *
                               in[nIds[r] * is + point];
                if (sum <= 0)
                    continue;
                par = p.kipar * n * pow(temperature, 2.5) /
                      (sqrt(A) * Z * Z * sum);
                perp  = p.kiperp * sqrt(A) * sum * n / (B2 * sqrt(temperature));
                cross = p.kicross * n * temperature / (Z * sqrt(B2));
            }
            T gradient[3] = {0, 0, 0};
            for (unsigned d = 0; d < dim; ++d)
            {
                const T gn =
                    electron ? gradNe[d] : q[(nIds[s] * dim + d) * qs + point];
                gradient[d] = ((T(2) / 3) * q[(energy * dim + d) * qs + point] -
                               temperature * gn) /
                              n;
            }
            for (unsigned d = 0; d < dim; ++d)
            {
                const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                T flux = cross * (b[j] * gradient[k] - b[k] * gradient[j]);
                for (unsigned r = 0; r < dim; ++r)
                    flux += ((par - perp) * b[d] * b[r] + (d == r ? perp : 0)) *
                            gradient[r];
                out[(energy * dim + d) * os + point] += scale * flux;
            }
        }
    }
};
} // namespace PENKNIFE::closures
