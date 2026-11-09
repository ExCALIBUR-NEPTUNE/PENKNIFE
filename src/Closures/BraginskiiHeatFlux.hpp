#pragma once

#include "Misc/Constants.hpp"
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <cmath>
#include <string>

namespace PENKNIFE::closures
{
// Heat-conduction policy; independent of the legacy Array-based closure.
template <typename TData> struct Braginskii
{
    static inline const std::string name = "Braginskii";

    static constexpr bool RequiresEnergy = true;

    struct Parameters
    {
        TData Nnorm, Tnorm, omega_c, mesh_length;
    };

    // Frequencies are in omega_c units, temperatures in Tnorm eV and densities
    // in Nnorm. No volume friction or energy-exchange sources are written here.
    NEK_HOSTDEVICE_INLINE static void Evaluate(
        const Parameters &p, const size_t point, const unsigned nDim,
        const unsigned nSpecies, const unsigned electronEnergyId,
        const unsigned *densityIds, const unsigned *energyIds,
        const TData *masses, const TData *charges, const size_t inStride,
        const size_t derivStride, const size_t outStride, const size_t bStride,
        const TData *in, const TData *deriv, const TData *b, TData *out,
        TData scale)
    {
        TData ne        = 0;
        TData gradNe[3] = {0, 0, 0};
        for (unsigned s = 0; s < nSpecies; ++s)
        {
            if (charges[s] == 0)
                continue;
            ne += charges[s] * in[densityIds[s] * inStride + point];
            for (unsigned d = 0; d < nDim; ++d)
                gradNe[d] +=
                    charges[s] *
                    deriv[(densityIds[s] * nDim + d) * derivStride + point];
        }
        // Padded Field lanes have zero density: do not evaluate log(0) there.
        if (ne == 0)
            return;
        const TData Te =
            (TData(2) / 3) * in[electronEnergyId * inStride + point] / ne;
        const TData collisionScale = TData(constants::c) /
                                     sqrt(TData(constants::m_p)) * p.Nnorm /
                                     TData(1e12) / p.omega_c;
        const TData conductionScale =
            scale * TData(constants::qeomp) * p.Tnorm /
            (p.omega_c * p.mesh_length * p.mesh_length);
        TData nuE = ElectronElectron(p, ne, Te) * collisionScale;

        for (unsigned s = 0; s < nSpecies; ++s)
        {
            if (charges[s] == 0)
                continue;
            const TData ni = in[densityIds[s] * inStride + point];
            if (ni == 0)
                continue;
            const TData Ti =
                (TData(2) / 3) * in[energyIds[s] * inStride + point] / ni;
            const TData nuEi =
                ElectronIon(p, ni, ne, Ti, Te, masses[s], charges[s]) *
                collisionScale;
            nuE += nuEi;
            TData nuI = TData(constants::m_e_m_p) * ne * nuEi / ni;
            // Sum ordered ion partners, including self collisions once.
            for (unsigned t = 0; t < nSpecies; ++t)
            {
                if (charges[t] == 0)
                    continue;
                const TData nt = in[densityIds[t] * inStride + point];
                if (nt == 0)
                    continue;
                const TData Tt =
                    (TData(2) / 3) * in[energyIds[t] * inStride + point] / nt;
                nuI += IonIon(p, ni, nt, Ti, Tt, masses[s], masses[t],
                              charges[s], charges[t]) *
                       collisionScale;
            }
            const TData kpar = nuI == 0 ? 0
                                        : conductionScale * TData(3.9) * ni *
                                              Ti / (masses[s] * nuI);
            AddHeatFlux(point, nDim, densityIds[s], energyIds[s], ni, Ti, kpar,
                        derivStride, outStride, bStride, deriv, b, out);
        }
        const TData kparE = nuE == 0 ? 0
                                     : conductionScale * TData(3.16) * ne * Te /
                                           (TData(constants::m_e_m_p) * nuE);
        TData gradTparE   = 0;
        for (unsigned d = 0; d < nDim; ++d)
        {
            const TData gradU =
                deriv[(electronEnergyId * nDim + d) * derivStride + point];
            const TData gradT = ((TData(2) / 3) * gradU - Te * gradNe[d]) / ne;
            gradTparE += b[d * bStride + point] * gradT;
        }
        for (unsigned d = 0; d < nDim; ++d)
            out[(electronEnergyId * nDim + d) * outStride + point] +=
                kparE * b[d * bStride + point] * gradTparE;
    }

private:
    NEK_HOSTDEVICE_INLINE static void AddHeatFlux(
        size_t point, unsigned nDim, unsigned nId, unsigned eId, TData n,
        TData T, TData kpar, size_t derivStride, size_t outStride,
        size_t bStride, const TData *deriv, const TData *b, TData *out)
    {
        TData gradTpar = 0;
        for (unsigned d = 0; d < nDim; ++d)
        {
            const TData gradU = deriv[(eId * nDim + d) * derivStride + point];
            const TData gradN = deriv[(nId * nDim + d) * derivStride + point];
            gradTpar += b[d * bStride + point] *
                        ((TData(2) / 3) * gradU - T * gradN) / n;
        }
        // Nektar's diffusion operator consumes +k grad(T), i.e. minus the
        // physical heat flux, to assemble +div(k grad(T)) in the RHS.
        for (unsigned d = 0; d < nDim; ++d)
            out[(eId * nDim + d) * outStride + point] +=
                kpar * b[d * bStride + point] * gradTpar;
    }

    NEK_HOSTDEVICE_INLINE static TData Denominator(TData speedSquared,
                                                   TData mass)
    {
        const TData pi          = TData(3.1415926535897932384626433832795);
        const TData x           = pi * speedSquared;
        const TData epsilonMass = TData(constants::epsilon_0) * mass;
        return TData(3) * x * sqrt(x) * epsilonMass * epsilonMass;
    }

    NEK_HOSTDEVICE_INLINE static TData ElectronElectron(const Parameters &p,
                                                        TData ne, TData Te)
    {
        const TData logTe = log(p.Tnorm * Te);
        const TData logLambda =
            TData(30.4) - TData(0.5) * log(ne * p.Nnorm) + TData(1.25) * logTe -
            sqrt(TData(1e-5) + (logTe - 2) * (logTe - 2) / 16);
        const TData v2 = 2 * p.Tnorm * Te / TData(constants::m_e_m_p);
        return 2 * ne * logLambda /
               Denominator(2 * v2, TData(constants::m_e_m_p));
    }

    NEK_HOSTDEVICE_INLINE static TData ElectronIon(const Parameters &p,
                                                   TData ni, TData ne, TData Ti,
                                                   TData Te, TData A, TData Z)
    {
        TData logLambda;
        if (Te * p.Tnorm < TData(0.1) || ni * p.Nnorm < TData(1e10) ||
            ne * p.Nnorm < TData(1e10))
            logLambda = 10;
        else if (Te < Ti * TData(constants::m_e_m_p) / A)
            logLambda = 23 - TData(0.5) * log(ni * p.Nnorm) +
                        TData(1.5) * log(Ti * p.Tnorm) - log(Z * Z * A);
        else if (Te * p.Tnorm < exp(TData(2)) * Z * Z)
            logLambda = 30 - TData(0.5) * log(ne * p.Nnorm) - log(Z) +
                        TData(1.5) * log(Te * p.Tnorm);
        else
            logLambda = 31 - TData(0.5) * log(ne * p.Nnorm) + log(Te * p.Tnorm);
        const TData ve2 = 2 * p.Tnorm * Te / TData(constants::m_e_m_p);
        const TData vi2 = 2 * p.Tnorm * Ti / A;
        // Preserve the electron-ion prefactor of the reference closure.
        return Z * Z * ni * logLambda * (1 + TData(constants::m_e_m_p)) /
               Denominator(ve2 + vi2, TData(constants::m_e_m_p));
    }

    NEK_HOSTDEVICE_INLINE static TData IonIon(const Parameters &p, TData ns,
                                              TData nt, TData Ts, TData Tt,
                                              TData As, TData At, TData Zs,
                                              TData Zt)
    {
        // Coulomb-log temperatures must be in eV, as in the thermal speeds.
        const TData TsEV = p.Tnorm * Ts, TtEV = p.Tnorm * Tt;
        const TData logLambda =
            TData(29.91) - TData(0.5) * log(p.Nnorm) -
            log(Zs * Zt * (As + At) / (As * TtEV + At * TsEV) *
                sqrt(ns * Zs * Zs / TsEV + nt * Zt * Zt / TtEV));
        const TData vs2 = 2 * TsEV / As, vt2 = 2 * TtEV / At;
        return Zs * Zs * Zt * Zt * nt * logLambda * (1 + As / At) /
               Denominator(vs2 + vt2, As);
    }
};
} // namespace PENKNIFE::closures
