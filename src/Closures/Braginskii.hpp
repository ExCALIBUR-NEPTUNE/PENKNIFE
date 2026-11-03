#ifndef BRAGINSKII_HPP
#define BRAGINSKII_HPP

#include "Closure.hpp"

namespace PENKNIFE
{

template <typename TData>
NEK_HOSTDEVICE_INLINE static TData CoulombLog_ii(TData Nnorm, TData ni1,
                                                 TData ni2, TData Ti1,
                                                 TData Ti2, TData A1, TData A2,
                                                 TData Z1, TData Z2)
{
    return 29.91 - log(sqrt(Nnorm)) -
           log((Z1 * Z2 * (A1 + A2)) / (A1 * Ti2 + A2 * Ti1) *
               sqrt(ni1 * Z1 * Z1 / Ti1 + ni2 * Z2 * Z2 / Ti2));
}
template <typename TData>
NEK_HOSTDEVICE_INLINE static TData CoulombLog_ee(TData Nnorm, TData Tnorm,
                                                 TData ne, TData Te)
{
    TData logTe = log(Tnorm * Te);
    return 30.4 - 0.5 * log(ne) - 0.5 * log(Nnorm) + (5. / 4) * logTe -
           sqrt(1e-5 + (logTe - 2) * (logTe - 2) / 16.);
}
template <typename TData>
NEK_HOSTDEVICE_INLINE static TData CoulombLog_ei(TData Nnorm, TData Tnorm,
                                                 TData ni, TData ne, TData Ti,
                                                 TData Te, TData Ai, TData Zi)
{
    if ((Te * Tnorm < 0.1) || (ni * Nnorm < 1e10) || (ne * Nnorm < 1e10))
        return 10;
    else if (Te < Ti * constants::m_e_m_p / Ai)
        return 23 - 0.5 * log(ni) + 1.5 * log(Ti) - log(Zi * Zi * Ai) -
               0.5 * log(Nnorm) + 1.5 * log(Tnorm);
    else if (Te * Tnorm < exp(2) * Zi * Zi)
        return 30.0 - 0.5 * log(ne) - log(Zi) + 1.5 * log(Te) -
               0.5 * log(Nnorm) + 1.5 * log(Tnorm);
    else
        return 31.0 - 0.5 * log(ne) + log(Te) - 0.5 * log(Nnorm) + log(Tnorm);
}

template <typename TData> constexpr TData BraginskiiCm(const TData Z)
{
    if (Z == 1)
        return 0.51;
    else if (Z == 2)
        return 0.44;
    else if (Z == 3)
        return 0.40;
    else
        return 0.38;
}

template <typename TData>
NEK_HOSTDEVICE_INLINE static TData PerpIonConductivity(TData n, TData T,
                                                       TData nu, TData m,
                                                       TData q, TData Bsq)
{
    TData Omegasq = constants::qeomp * constants::qeomp * Bsq;
    return 2.0 * n * T * nu / (m * Omegasq);
}

template <typename TData>
NEK_HOSTDEVICE_INLINE static TData PerpElectronConductivity(TData n, TData T,
                                                            TData nu, TData Bsq)
{
    TData Omegasq = (constants::e / constants::m_e_si) *
                    (constants::e / constants::m_e_si) * Bsq;

    return (sqrt(2.0) + 3.25) * n * T * nu / (constants::m_e_m_p * Omegasq);
}

template <typename TData>
NEK_HOSTDEVICE_INLINE static TData CrossIonConductivity(TData n, TData T,
                                                        TData q, TData B)
{
    return 2.5 * n * T / (q * B);
}

template <typename TData>
NEK_HOSTDEVICE_INLINE static TData CrossElectronConductivity(TData n, TData T,
                                                             TData B)
{
    return 2.5 * n * T / B;
}

template <typename TData>
NEK_HOSTDEVICE_INLINE static void DoBraginskiiClosure(
    const unsigned i, const unsigned n_species, const unsigned int *n_ids,
    const unsigned int *v_ids, const unsigned int *e_ids, const TData *mass,
    const TData *charge, const TData Nnorm, const TData TNorm,
    const size_t inVecStride, const size_t freqStride,
    const size_t frictionStride, const TData *in, TData *frequencies,
    TData *frictions)
{
    const size_t ee_off = ee_idx * inVecStride + i;
    const TData Te      = (2. / 3) * in[ee_off] / ne[i];

    const TData v1sq  = 2 * TNorm * Te / constants::m_e_m_p;
    TData coulomb_log = CoulombLog_ee(Nnorm, Tnorm, ne[p], Te);

    // Electon collision frequency
    TData nu = ne[p] * coulomb_log * 2 /
               (3 * pow(M_PI * 2 * v1sq, 1.5) *
                pow(constants::epsilon_0 * constants::m_e_m_p, 2));
    nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
    // nu in s^-1
    nu /= omega_c;
    // this->nu_ee[p] = nu;
    //  nu in gyrofrequencies

    // this->nu_e[p] = this->nu_ee[p];
    frequencies[n_species * freqStride + i] = nu;

    for (unsigned int s = 0; s < n_species; ++s)
    {
        TData Z                   = charge[s];
        TData A                   = mass[s];
        const unsigned int ni_idx = n_ids[s];
        const unsigned int ei_idx = v_ids[s];

        const size_t n_off = ni_idx * inVecStride + i;
        const size_t e_off = ei_idx * inVecStride + i;

        const TData Ti   = (2. / 3) * in[e_off] / in[n_off];
        const TData Te   = (2. / 3) * in[ee_off] / ne[i];
        const TData vesq = 2 * TNorm * Te / constants::m_e_m_p;
        const TData visq = 2 * TNorm * Ti / A;
        TData coulomb_log =
            CoulombLog_ei(Nnorm, Tnorm, in[n_off], ne[i], Ti, Te, A, Z);
        // Collision frequency
        TData nu = Z * Z * in[n_off] * coulomb_log * (1. + constants::m_e_m_p) /
                   (3 * pow(M_PI * (vesq + visq), 1.5) *
                    pow(constants::epsilon_0 * constants::m_e_m_p, 2));
        nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
        nu /= omega_c;

        // this->nu_ei[s][p] = nu;
        //  this->nu_e[p] += nu;
        frequencies[n_species * freqStride + i] += nu;
        // this->nu_i[s][p] = constants::m_e_m_p * ne[i] * nu / in[n_off];
        frequencies[s * freqStride + i] =
            constants::m_e_m_p * ne[i] * nu / in[n_off];

        for (unsigned int s2 = 0; s2 < n_species; ++s2)
        {
            if (s2 > s)
                break;
            TData Z2                   = charge[s2];
            TData A2                   = mass[s2];
            const unsigned int ni_idx2 = n_ids[s2];
            const unsigned int ei_idx2 = v_ids[s2];
            const size_t n_off2        = ni_idx2 * inVecStride + i;
            const size_t e_off2        = ei_idx2 * inVecStride + i;

            const TData Ti2 = (2. / 3) * in[e_off2] / in[n_off2];

            TData coulomb_log = CoulombLog_ii(Nnorm, in[n_off], in[n_off2], Ti,
                                              Ti2, A, A2, Z, Z2);

            const TData v1sq = 2 * Tnorm * Ti / A;
            const TData v2sq = 2 * Tnorm * Ti2 / A2;
            TData nu         = Z * Z * Z2 * Z2 * in[n_off2] * coulomb_log *
                               (1. + A / A2) /
                               (3 * pow(M_PI * (v1sq + v2sq), 1.5) *
                                pow(A * constants::epsilon_0, 2));
            nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
            nu /= omega_c;

            // this->nu_ii[std::make_pair(s, s2)][p] = nu;
            // this->nu_i[s][p] += nu;
            frequencies[s * freqStride + i] += nu;
            // this->nu_i[s2][p] += (A / A2) * nu * in[n_off] / in[n_off2];
            frequencies[s2 * freqStride + i] +=
                (A / A2) * nu * in[n_off] / in[n_off2];
        }
    }
}

/**
 *
 */
class Braginskii : public Closure
{
public:
    friend class MemoryManager<Braginskii>;

    /// Creates an instance of this class
    static ClosureSharedPtr create(const std::weak_ptr<PlasmaSystem> &pSystem,
                                   const int spaceDim)
    {
        ClosureSharedPtr p =
            MemoryManager<Braginskii>::AllocateSharedPtr(pSystem, spaceDim);
        return p;
    }

    /// Name of the class
    static std::string className;

private:
    void v_EvaluateClosure(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        const Array<OneD, NekDouble> &ne,
        const Array<OneD, NekDouble> &ve) override;

    void CalcCollisionFrequencies(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, NekDouble> &ne);

    Braginskii(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim);

    ~Braginskii() = default;

    LibUtilities::Field<TData, FieldState::Phys> frequencies;

    NekDouble k_ci;
    NekDouble k_ce;
};

} // namespace PENKNIFE

#endif