#include "Braginskii.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

namespace PENKNIFE
{
std::string Braginskii::className = GetClosureFactory().RegisterCreatorFunction(
    "Braginskii", Braginskii::create, "Braginskii closure system");

Braginskii::Braginskii(const std::weak_ptr<PlasmaSystem> &pSystem,
                       const int spaceDim)
    : Closure(pSystem, spaceDim)
{
    k_ci        = 3.9;
    k_ce        = 3.16;
    this->nu_ee = Array<OneD, NekDouble>(this->n_pts);
    this->nu_e  = Array<OneD, NekDouble>(this->n_pts);
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        this->nu_i[s]  = Array<OneD, NekDouble>(this->n_pts);
        this->nu_ei[s] = Array<OneD, NekDouble>(this->n_pts);
        for (const auto &[s2, v2] : m_system.lock()->GetIons())
        {
            this->nu_ii[std::make_pair(s, s2)] =
                Array<OneD, NekDouble>(this->n_pts);
        }
    }
}

inline double CoulombLog_ii(double Nnorm, double ni1, double ni2, double Ti1,
                            double Ti2, double A1, double A2, double Z1,
                            double Z2)
{
    return 29.91 - log(sqrt(Nnorm)) -
           log((Z1 * Z2 * (A1 + A2)) / (A1 * Ti2 + A2 * Ti1) *
               sqrt(ni1 * Z1 * Z1 / Ti1 + ni2 * Z2 * Z2 / Ti2));
}

inline double CoulombLog_ee(double Nnorm, double Tnorm, double ne, double Te)
{
    double logTe = log(Tnorm * Te);
    return 30.4 - 0.5 * log(ne) - 0.5 * log(Nnorm) + (5. / 4) * logTe -
           sqrt(1e-5 + (logTe - 2) * (logTe - 2) / 16.);
}

inline double CoulombLog_ei(double Nnorm, double Tnorm, double ni, double ne,
                            double Ti, double Te, double Ai, double Zi)
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

void Braginskii::CalcCollisionFrequencies(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, NekDouble> &ne)
{
    for (int p = 0; p < this->n_pts; ++p)
    {
        const double v1sq =
            2 * in_arr[ee_idx][p] * constants::e / constants::m_e_si;
        double coulomb_log =
            CoulombLog_ee(Nnorm, Tnorm, ne[p], in_arr[ee_idx][p]);

        // Electon collision frequency
        double nu = pow(constants::e, 4) * ne[p] * coulomb_log * 2 /
                    (3 * pow(M_PI * 2 * v1sq, 1.5) *
                     pow(constants::epsilon_0_si * constants::m_e_si, 2));
        this->nu_ee[p] = nu / omega_c;
        this->nu_e[p]  = this->nu_ee[p];
    }
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        double Z   = v.charge;
        double A   = v.mass;
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int ei_idx = v.fields.at(field_to_index.at("e"));
        for (int p = 0; p < this->n_pts; ++p)
        {
            const double vesq =
                2 * in_arr[ei_idx][p] * constants::e / constants::m_e_si;
            const double visq =
                2 * in_arr[ei_idx][p] * constants::e / (A * constants::m_p_si);
            double coulomb_log =
                CoulombLog_ei(Nnorm, Tnorm, in_arr[ni_idx][p], ne[p],
                              in_arr[ei_idx][p], in_arr[ee_idx][p], A, Z);
            // Collision frequency
            double nu = Z * Z * pow(constants::e, 4) * in_arr[ni_idx][p] *
                        coulomb_log * (1. + constants::m_e_m_p) /
                        (3 * pow(M_PI * (vesq + visq), 1.5) *
                         pow(constants::epsilon_0_si * constants::m_e_si, 2));
            nu /= omega_c;
            this->nu_ei[s][p] = nu;
            this->nu_e[p] += nu;
            this->nu_i[s][p] =
                constants::m_e_m_p * ne[p] * nu / in_arr[ni_idx][p];
        }
    }
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        double Z   = v.charge;
        double A   = v.mass;
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int ei_idx = v.fields.at(field_to_index.at("e"));
        for (const auto &[s2, v2] : m_system.lock()->GetIons())
        {
            if (s2 > s)
                break;
            double Z2   = v2.charge;
            double A2   = v2.mass;
            int ni_idx2 = v2.fields.at(field_to_index.at("n"));
            int ei_idx2 = v2.fields.at(field_to_index.at("e"));
            for (int p = 0; p < this->n_pts; ++p)
            {
                double coulomb_log = CoulombLog_ii(
                    Nnorm, in_arr[ni_idx][p], in_arr[ni_idx2][p],
                    in_arr[ei_idx][p], in_arr[ei_idx2][p], A, A2, Z, Z2);

                const double v1sq = 2 * in_arr[ei_idx][p] * constants::e /
                                    (A * constants::m_p_si);
                const double v2sq = 2 * in_arr[ei_idx2][p] * constants::e /
                                    (A2 * constants::m_p_si);
                double nu =
                    Z * Z * Z2 * Z2 * pow(constants::e, 4) *
                    in_arr[ni_idx2][p] * coulomb_log * (1. + A / A2) /
                    (3 * pow(M_PI * (v1sq + v2sq), 1.5) *
                     pow(constants::epsilon_0_si * A * constants::m_p_si, 2));
                nu /= omega_c;
                this->nu_ii[std::make_pair(s, s2)][p] = nu;
                this->nu_i[s][p] += nu;
                this->nu_i[s2][p] +=
                    (A / A2) * nu * in_arr[ni_idx][p] / in_arr[ni_idx][p];
            }
        }
    }
}

constexpr inline double BraginskiiCm(const double Z)
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

inline double PerpIonConductivity(double n, double T, double nu, double m,
                                  double q, double Bsq)
{
    return 2 * n * T * m * nu / (q * q * Bsq);
}

inline double PerpElectronConductivity(double n, double T, double nu,
                                       double Bsq)
{
    return (sqrt(2.0) + 3.25) * n * T * constants::m_e_m_p * nu / Bsq;
}

inline double CrossIonConductivity(double n, double T, double q, double B)
{
    return 2.5 * n * T / (q * B);
}

inline double CrossElectronConductivity(double n, double T, double B)
{
    return 2.5 * n * T / B;
}

void Braginskii::v_EvaluateClosure(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
    Array<OneD, Array<OneD, NekDouble>> &frictions,
    const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve)
{
    CalcCollisionFrequencies(values, ne);

    // Heat fluxes
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int ei_idx = v.fields.at(field_to_index.at("e"));

        for (int p = 0; p < this->n_pts; ++p)
        {
            double kpar =
                k_ci * values[ei_idx][p] * values[ni_idx][p] / nu_i[s][p];
            double kperp =
                PerpIonConductivity(values[ni_idx][p], values[ei_idx][p],
                                    nu_i[s][p], v.mass, v.charge, mag_B[p]);
            double kcross = CrossIonConductivity(
                values[ni_idx][p], values[ei_idx][p], v.charge, sqrt(mag_B[p]));

            for (unsigned int i = 0; i < m_spacedim; ++i)
            {
                for (unsigned int j = 0; j < m_spacedim; ++j)
                {
                    double D = (kpar - kperp) * b_unit[i][p] * b_unit[j][p];
                    if (i == j)
                    {
                        D += kperp;
                    }
                    fluxes[i][ei_idx][p] += D * grads[j][ei_idx][p];
                }
            }
            if (m_spacedim == 3)
            {
                fluxes[0][ei_idx][p] +=
                    kcross * (b_unit[1][p] * grads[2][ei_idx][p] -
                              b_unit[2][p] * grads[1][ei_idx][p]);
                fluxes[1][ei_idx][p] +=
                    kcross * (b_unit[2][p] * grads[0][ei_idx][p] -
                              b_unit[0][p] * grads[2][ei_idx][p]);
                fluxes[2][ei_idx][p] +=
                    kcross * (b_unit[0][p] * grads[1][ei_idx][p] -
                              b_unit[1][p] * grads[0][ei_idx][p]);
            }
            else
            {
                fluxes[0][ei_idx][p] +=
                    -kcross * b_unit[2][p] * grads[1][ei_idx][p];
                fluxes[1][ei_idx][p] +=
                    kcross * b_unit[2][p] * grads[0][ei_idx][p];
            }
        }
    }

    for (int p = 0; p < this->n_pts; ++p)
    {
        double kpar  = k_ce * values[ee_idx][p] * ne[p] / nu_e[p];
        double kperp = PerpElectronConductivity(ne[p], values[ee_idx][p],
                                                nu_e[p], mag_B[p]);
        double kcross =
            CrossElectronConductivity(ne[p], values[ee_idx][p], sqrt(mag_B[p]));

        for (unsigned int i = 0; i < m_spacedim; ++i)
        {
            for (unsigned int j = 0; j < m_spacedim; ++j)
            {
                double D = (kpar - kperp) * b_unit[i][p] * b_unit[j][p];
                if (i == j)
                {
                    D += kperp;
                }
                fluxes[i][ee_idx][p] += D * grads[j][ee_idx][p];
            }
        }
        if (m_spacedim == 3)
        {
            fluxes[0][ee_idx][p] +=
                kcross * (b_unit[1][p] * grads[2][ee_idx][p] -
                          b_unit[2][p] * grads[1][ee_idx][p]);
            fluxes[1][ee_idx][p] +=
                kcross * (b_unit[2][p] * grads[0][ee_idx][p] -
                          b_unit[0][p] * grads[2][ee_idx][p]);
            fluxes[2][ee_idx][p] +=
                kcross * (b_unit[0][p] * grads[1][ee_idx][p] -
                          b_unit[1][p] * grads[0][ee_idx][p]);
        }
        else
        {
            fluxes[0][ee_idx][p] +=
                -kcross * b_unit[2][p] * grads[1][ee_idx][p];
            fluxes[1][ee_idx][p] += kcross * b_unit[2][p] * grads[0][ee_idx][p];
        }
    }

    // Friction and heat exchange
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int vi_idx = v.fields.at(field_to_index.at("v"));
        int ei_idx = v.fields.at(field_to_index.at("e"));
        double mu  = constants::m_e_m_p / (v.mass + constants::m_e_m_p);
        double cm  = BraginskiiCm(v.charge) * v.mass;
        for (int p = 0; p < this->n_pts; ++p)
        {
            // Friction momentum (electron to ion)
            double Fie = cm * nu_ei[s][p] * (ve[p] - values[vi_idx][p]);
            frictions[vi_idx][p] = Fie;

            // Frictional heating (electron to ion)
            double Qie = mu * Fie * (ve[p] - values[vi_idx][p]);

            // Heat exchange (electron to ion)
            Qie += mu * nu_ei[s][p] * values[ni_idx][p] *
                   (values[ee_idx][p] - values[ei_idx][p]);
            frictions[ei_idx][p] = Qie;
            frictions[ee_idx][p] = -Qie;
        }

        for (const auto &[s2, v2] : m_system.lock()->GetIons())
        {
            if (s2 >= s)
                break;
            int ni_idx2 = v2.fields.at(field_to_index.at("n"));
            int vi_idx2 = v2.fields.at(field_to_index.at("v"));
            int ei_idx2 = v2.fields.at(field_to_index.at("e"));
            double mu   = v2.mass / (v.mass + v2.mass);

            for (int p = 0; p < this->n_pts; ++p)
            {
                // Friction momentum (ion2 to ion)
                double Fii = 1.0 * nu_ii[std::make_pair(s, s2)][p] * v.mass *
                             (values[vi_idx2][p] - values[vi_idx][p]);
                frictions[vi_idx][p] += Fii;
                frictions[vi_idx2][p] -= Fii;

                // Frictional heating (ion2 to ion)
                double Qii =
                    mu * Fii * (values[vi_idx2][p] - values[vi_idx][p]);

                // Heat exchange (ion2 to ion)
                Qii += mu * nu_ii[std::make_pair(s, s2)][p] *
                       values[ni_idx][p] *
                       (values[ei_idx2][p] - values[ei_idx][p]);

                frictions[ei_idx][p] += Qii;
                frictions[ei_idx2][p] -= Qii;
            }
        }
    }
}
} // namespace PENKNIFE
