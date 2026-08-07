#include "Braginskii.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

namespace PENKNIFE
{
std::string Braginskii::className = GetClosureFactory().RegisterCreatorFunction(
    "Braginskii", Braginskii::create, "Braginskii closure system");

/**
 * @brief Constructor for the Braginskii closure.
 */
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

/**
 * @brief Constructor for the Braginskii closure.
 * @param values Physical field values
 * @param ne Electron density
 */
void Braginskii::v_CollisionFrequencies(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, NekDouble> &ne)
{
    for (int p = 0; p < this->n_pts; ++p)
    {
        double Te          = (2. / 3.) * values[ee_idx][p] / ne[p];
        const double v1sq  = 2 * Te / constants::m_e_m_p;
        double coulomb_log = CoulombLog_ee(Nnorm, Tnorm, ne[p], Te);

        // Electon collision frequency
        double nu = ne[p] * coulomb_log * 2 /
                    (3 * pow(M_PI * 2 * v1sq, 1.5) *
                     pow(constants::epsilon_0 * constants::m_e_m_p, 2));
        nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
        // nu in s^-1
        nu /= omega_c;
        this->nu_ee[p] = nu;
        // nu in gyrofrequencies

        this->nu_e[p] = this->nu_ee[p];
    }
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        double Z   = v.charge;
        double A   = v.mass;
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int ei_idx = v.fields.at(field_to_index.at("e"));
        for (int p = 0; p < this->n_pts; ++p)
        {
            double Ti = (2. / 3.) * values[ei_idx][p] / values[ni_idx][p];
            double Te = (2. / 3.) * values[ee_idx][p] / ne[p];

            const double vesq  = 2 * Ti / constants::m_e_m_p;
            const double visq  = 2 * Ti / A;
            double coulomb_log = CoulombLog_ei(Nnorm, Tnorm, values[ni_idx][p],
                                               ne[p], Ti, Te, A, Z);
            // Collision frequency
            double nu = Z * Z * values[ni_idx][p] * coulomb_log *
                        (1. + constants::m_e_m_p) /
                        (3 * pow(M_PI * (vesq + visq), 1.5) *
                         pow(constants::epsilon_0 * constants::m_e_m_p, 2));
            nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
            nu /= omega_c;

            this->nu_ei[s][p] = nu;
            this->nu_e[p] += nu;
            this->nu_i[s][p] =
                constants::m_e_m_p * ne[p] * nu / values[ni_idx][p];


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
                double Ti = (2. / 3.) * values[ei_idx][p] / values[ni_idx][p];
                double Ti2 =
                    (2. / 3.) * values[ei_idx2][p] / values[ni_idx2][p];

                double coulomb_log =
                    CoulombLog_ii(Nnorm, values[ni_idx][p], values[ni_idx2][p],
                                  Ti, Ti2, A, A2, Z, Z2);

                const double v1sq = 2 * Tnorm * Ti / A;
                const double v2sq = 2 * Tnorm * Ti2 / A2;
                double nu = Z * Z * Z2 * Z2 * values[ni_idx2][p] * coulomb_log *
                            (1. + A / A2) /
                            (3 * pow(M_PI * (v1sq + v2sq), 1.5) *
                             pow(A * constants::epsilon_0, 2));
                nu *= (constants::c / (sqrt(constants::m_p))) * Nnorm / 1e12;
                nu /= omega_c;

                this->nu_ii[std::make_pair(s, s2)][p] = nu;
                this->nu_i[s][p] += nu;
                this->nu_i[s2][p] +=
                    (A / A2) * nu * values[ni_idx][p] / values[ni_idx][p];
            }
        }
    }
}

void Braginskii::v_EvaluateConductivity(const Array<OneD, NekDouble> &ne,
                                        Array<OneD, NekDouble> &sigma)
{
    
    for (int p = 0; p < this->n_pts; ++p)
    {

        sigma[p] = 1.96 * ne[p] / (Bnorm * this->nu_e[p] * constants::m_e_m_p);

        // std::cout<<"sigma "<<sigma[p]<<"\n";
    }
    //std::cout<<"sigma "<<sigma[1223]<<"nu "<<nu_e[1223]<<"\n";

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
    double Omegasq = constants::qeomp * constants::qeomp * Bsq;
    return 2.0 * n * T * nu / (m * Omegasq);
}

inline double PerpElectronConductivity(double n, double T, double nu,
                                       double Bsq)
{
    double Omegasq = (constants::e / constants::m_e_si) *
                     (constants::e / constants::m_e_si) * Bsq;

    return (sqrt(2.0) + 3.25) * n * T * nu / (constants::m_e_m_p * Omegasq);
}

inline double CrossIonConductivity(double n, double T, double q, double B)
{
    return 2.5 * n * T / (q * B);
}

inline double CrossElectronConductivity(double n, double T, double B)
{
    return 2.5 * n * T / B;
}

/**
 * @brief Evaluate the Braginskii closure.
 * @param values Physical field values
 * @param grads Physical field gradients
 * @param[out] fluxes Heat fluxes
 * @param ne Electron density
 */
void Braginskii::v_EvaluateHeatFlux(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    int ti_idx = 0;
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int ei_idx = v.fields.at(field_to_index.at("e"));

        for (int p = 0; p < this->n_pts; ++p)
        {
            double Ti = (2. / 3.) * values[ei_idx][p] / values[ni_idx][p];

            double kpar = scaling * (constants::qeomp / omega_c) * k_ci *
                          Tnorm * Ti * values[ni_idx][p] /
                          (v.mass * nu_i[s][p]);

            double kperp =
                scaling * constants::qeomp * omega_c * Tnorm *
                PerpIonConductivity(values[ni_idx][p], Ti, nu_i[s][p], v.mass,
                                    v.charge, mag_B[p]);
            double kcross = scaling * Tnorm *
                            CrossIonConductivity(values[ni_idx][p], Ti,
                                                 v.charge, sqrt(mag_B[p]));

            // std::cout << "ion kpar = " << kpar << " ion kperp = " << kperp
            //           << " ion kcross = " << kcross << "\n";

            for (unsigned int i = 0; i < m_spacedim; ++i)
            {
                for (unsigned int j = 0; j < m_spacedim; ++j)
                {
                    double D = (kpar - kperp) * b_unit[i][p] * b_unit[j][p];
                    if (i == j)
                    {
                        D += kperp;
                    }
                    fluxes[i][ei_idx][p] += D * grads[j][ti_idx][p];
                }
            }
            if (m_spacedim == 3)
            {
                fluxes[0][ei_idx][p] +=
                    kcross * (b_unit[1][p] * grads[2][ti_idx][p] -
                              b_unit[2][p] * grads[1][ti_idx][p]);
                fluxes[1][ei_idx][p] +=
                    kcross * (b_unit[2][p] * grads[0][ti_idx][p] -
                              b_unit[0][p] * grads[2][ti_idx][p]);
                fluxes[2][ei_idx][p] +=
                    kcross * (b_unit[0][p] * grads[1][ti_idx][p] -
                              b_unit[1][p] * grads[0][ti_idx][p]);
            }
            else
            {
                fluxes[0][ei_idx][p] +=
                    -kcross * b_unit[2][p] * grads[1][ti_idx][p];
                fluxes[1][ei_idx][p] +=
                    kcross * b_unit[2][p] * grads[0][ti_idx][p];
            }
        }
        ti_idx++;
    }

    int te_idx = ti_idx;
    for (int p = 0; p < this->n_pts; ++p)
    {
        double Te = (2. / 3.) * values[ee_idx][p] / ne[p];

        double kpar   = scaling * (constants::qeomp / omega_c) * k_ce * Tnorm *
                        Te * ne[p] / (constants::m_e_m_p * nu_e[p]);
        double kperp  = scaling * constants::qeomp * omega_c * Tnorm *
                        PerpElectronConductivity(ne[p], Te, nu_e[p], mag_B[p]);
        double kcross = scaling * Tnorm *
                        CrossElectronConductivity(ne[p], Te, sqrt(mag_B[p]));

        // std::cout << "e kpar = " << kpar << " e kperp = " << kperp
        //           << " e kcross = " << kcross << "\n";

        for (unsigned int i = 0; i < m_spacedim; ++i)
        {
            for (unsigned int j = 0; j < m_spacedim; ++j)
            {
                double D = (kpar - kperp) * b_unit[i][p] * b_unit[j][p];
                if (i == j)
                {
                    D += kperp;
                }
                fluxes[i][ee_idx][p] += D * grads[j][te_idx][p];
            }
        }
        if (m_spacedim == 3)
        {
            fluxes[0][ee_idx][p] +=
                kcross * (b_unit[1][p] * grads[2][te_idx][p] -
                          b_unit[2][p] * grads[1][te_idx][p]);
            fluxes[1][ee_idx][p] +=
                kcross * (b_unit[2][p] * grads[0][te_idx][p] -
                          b_unit[0][p] * grads[2][te_idx][p]);
            fluxes[2][ee_idx][p] +=
                kcross * (b_unit[0][p] * grads[1][te_idx][p] -
                          b_unit[1][p] * grads[0][te_idx][p]);
        }
        else
        {
            fluxes[0][ee_idx][p] +=
                -kcross * b_unit[2][p] * grads[1][te_idx][p];
            fluxes[1][ee_idx][p] += kcross * b_unit[2][p] * grads[0][te_idx][p];
        }
    }
}

void Braginskii::v_EvaluateThermalForce(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, NekDouble>> &force)
{
    int ti_idx = 0;
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        for (int p = 0; p < this->n_pts; ++p)
        {
            for (unsigned int d = 0; d < m_spacedim; ++d)
            {
                force[ti_idx][p] -= 0.71 * ne[p] * v.charge * v.charge *
                                    b_unit[d][p] * grads[d][ti_idx][p];
            }
        }
        ti_idx++;
    }
}

void Braginskii::v_EvaluateFrictionHeating(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ge,
    Array<OneD, Array<OneD, NekDouble>> &frictions,
    Array<OneD, Array<OneD, NekDouble>> &heats)
{
    // Friction and heat exchange
    int te_idx = frictions.size() - 1;
    int ti_idx = 0;
    for (const auto &[s, v] : m_system.lock()->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index.at("n"));
        int vi_idx = v.fields.at(field_to_index.at("v"));
        int ei_idx = v.fields.at(field_to_index.at("e"));
        double mu  = 1. / (v.mass + constants::m_e_m_p);
        double cm  = BraginskiiCm(v.charge) * v.mass;

        for (int p = 0; p < this->n_pts; ++p)
        {
            double vi = values[vi_idx][p] / (values[ni_idx][p] * v.mass);
            double ve = ge[p];
            double Ti = (2. / 3.) * values[ei_idx][p] / values[ni_idx][p];
            double Te = (2. / 3.) * values[ee_idx][p] / ne[p];
            // Friction momentum (ion to electron)
            double Fei = cm * omega_c * nu_ei[s][p] * constants::m_e_m_p *
                         ne[p] * (vi - ve);
            frictions[ti_idx][p] -= Fei;
            frictions[te_idx][p] += Fei;

            // Frictional heating (ion to electron)
            double Wei = mu * Fei * (vi - ve);

            // convert to Tnorm eV
            Wei *= mesh_length * mesh_length / (Tnorm * constants::qeomp);
            heats[ti_idx][p] += constants::m_e_m_p * Wei;
            heats[te_idx][p] += v.mass * Wei;

            // Heat exchange (ion to electron)
            double Qei = 3 * mu * omega_c * nu_ei[s][p] * constants::m_e_m_p *
                         ne[p] * (Ti - Te);
            heats[ti_idx][p] -= Qei;
            heats[te_idx][p] += Qei;
            // std::cout << "nu_ei " << omega_c * nu_ei[s][p] << " vi friction "
            //           << frictions[ti_idx][p] << " ei heat " <<
            //           heats[ti_idx][p]
            //           << " e friction " << frictions[te_idx][p] << " Qei "
            //           << Qei << " Fei " << Fei << " Wei " << Wei << " ve " <<
            //           ve
            //           << " vi " << vi << " Te " << Te << " Ti " << Ti << "\n
            //           ";
        }

        int ti_idx2 = 0;
        for (const auto &[s2, v2] : m_system.lock()->GetIons())
        {
            if (ti_idx2 >= ti_idx)
                break;
            int ni_idx2 = v2.fields.at(field_to_index.at("n"));
            int vi_idx2 = v2.fields.at(field_to_index.at("v"));
            int ei_idx2 = v2.fields.at(field_to_index.at("e"));
            double mu   = 1. / (v.mass + v2.mass);

            for (int p = 0; p < this->n_pts; ++p)
            {
                double vi = values[vi_idx][p] / (values[ni_idx][p] * v.mass);
                double vi2 =
                    values[vi_idx2][p] / (values[ni_idx2][p] * v2.mass);

                double Ti = (2. / 3.) * values[ei_idx][p] / values[ni_idx][p];
                double Ti2 =
                    (2. / 3.) * values[ei_idx2][p] / values[ni_idx2][p];

                // Friction momentum (ion2 to ion)
                double Fii = 1.0 * omega_c * nu_ii[std::make_pair(s, s2)][p] *
                             v.mass * values[ni_idx][p] * (vi2 - vi);
                frictions[ti_idx][p] += Fii;
                frictions[ti_idx2][p] -= Fii;

                // Frictional heating (ion2 to ion)
                double Wii = mu * Fii * (vi2 - vi);

                // convert to Tnorm eV
                Wii *= mesh_length * mesh_length / (Tnorm * constants::qeomp);
                heats[ti_idx][p] += v2.mass * Wii;
                heats[ti_idx2][p] += v.mass * Wii;

                // Heat exchange (ion2 to ion)
                double Qii = 3 * mu * omega_c *
                             nu_ii[std::make_pair(s, s2)][p] * v.mass *
                             values[ni_idx][p] * (Ti2 - Ti);

                heats[ti_idx][p] += Qii;
                heats[ti_idx2][p] -= Qii;
            }
            ti_idx2++;
        }
        ti_idx++;
    }
}
} // namespace PENKNIFE
