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
    for (const auto s : m_system.lock()->GetIons())
    {
        this->nu_i[s]  = Array<OneD, NekDouble>(this->n_pts);
        this->nu_ei[s] = Array<OneD, NekDouble>(this->n_pts);
        for (const auto s2 : m_system.lock()->GetIons())
        {
            this->nu_ii[std::make_pair(s, s2)] =
                Array<OneD, NekDouble>(this->n_pts);
        }
    }
}



/**
 * @brief Evaluate the Braginskii closure.
 * @param values Physical field values
 * @param grads Physical field gradients
 * @param[out] fluxes Heat fluxes
 * @param[out] frictions Collisional frictions
 * @param ne Electron density
 * @param ve Electron velocity
 */
void Braginskii::v_EvaluateClosure(
    const Array<OneD, Array<OneD, NekDouble>> &values,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
    Array<OneD, Array<OneD, NekDouble>> &frictions,
    const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve)
{
    CalcCollisionFrequencies(values, ne);

    // Heat fluxes
    for (const auto s : m_system.lock()->GetIons())
    {
        const double mass_s   = m_system.lock()->GetMasses().at(s);
        const double charge_s = m_system.lock()->GetCharges().at(s);
        int ni_idx            = m_system.lock()->GetDensityId(s);
        int ei_idx            = m_system.lock()->GetEnergyId(s);

        for (int p = 0; p < this->n_pts; ++p)
        {
            double kpar = scaling * (constants::qeomp / omega_c) * k_ci *
                          Tnorm * values[ei_idx][p] * values[ni_idx][p] /
                          (mass_s * nu_i[s][p]);

            double kperp =
                scaling * constants::qeomp * omega_c * Tnorm *
                PerpIonConductivity(values[ni_idx][p], values[ei_idx][p],
                                    nu_i[s][p], mass_s, charge_s, mag_B[p]);
            double kcross =
                scaling * Tnorm *
                CrossIonConductivity(values[ni_idx][p], values[ei_idx][p],
                                     charge_s, sqrt(mag_B[p]));

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
        double kpar  = scaling * (constants::qeomp / omega_c) * k_ce * Tnorm *
                       values[ee_idx][p] * ne[p] /
                       (constants::m_e_m_p * nu_e[p]);
        double kperp = scaling * constants::qeomp * omega_c * Tnorm *
                       PerpElectronConductivity(ne[p], values[ee_idx][p],
                                                nu_e[p], mag_B[p]);
        double kcross =
            scaling * Tnorm *
            CrossElectronConductivity(ne[p], values[ee_idx][p], sqrt(mag_B[p]));

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
    for (const auto s : m_system.lock()->GetIons())
    {
        const double mass_s   = m_system.lock()->GetMasses().at(s);
        const double charge_s = m_system.lock()->GetCharges().at(s);
        int ni_idx            = m_system.lock()->GetDensityId(s);
        int vi_idx            = m_system.lock()->GetMomentumId(s);
        int ei_idx            = m_system.lock()->GetEnergyId(s);
        double mu             = 1. / (mass_s + constants::m_e_m_p);
        double cm             = BraginskiiCm(charge_s) * mass_s;
        for (int p = 0; p < this->n_pts; ++p)
        {
            // Friction momentum (electron to ion)
            double Fie = cm * omega_c * nu_ei[s][p] * values[ni_idx][p] *
                         (ve[p] / (ne[p] * constants::m_e_m_p) -
                          values[vi_idx][p] / (values[ni_idx][p] * mass_s));
            frictions[vi_idx][p] = Fie;
            // No electron momentum equation

            // Frictional heating (electron to ion)
            double Wie = mu * Fie *
                         (ve[p] / (ne[p] * constants::m_e_m_p) -
                          values[vi_idx][p] / (values[ni_idx][p] * mass_s));

            // convert to Tnorm eV
            Wie *= mesh_length * mesh_length / (Tnorm * constants::qeomp);
            frictions[ei_idx][p] += constants::m_e_m_p * Wie;
            frictions[ee_idx][p] += mass_s * Wie;

            // Heat exchange (electron to ion)
            double Qie = 3 * mu * omega_c * nu_ei[s][p] * mass_s *
                         values[ni_idx][p] *
                         (values[ee_idx][p] - values[ei_idx][p]);
            frictions[ei_idx][p] += Qie;
            frictions[ee_idx][p] -= Qie;
        }

        for (const auto s2 : m_system.lock()->GetIons())
        {
            const double mass_s2 = m_system.lock()->GetMasses().at(s2);
            if (s2 >= s)
                break;
            int ni_idx2 = m_system.lock()->GetDensityId(s2);
            int vi_idx2 = m_system.lock()->GetMomentumId(s2);
            int ei_idx2 = m_system.lock()->GetEnergyId(s2);
            double mu   = 1. / (mass_s + mass_s2);

            for (int p = 0; p < this->n_pts; ++p)
            {
                // Friction momentum (ion2 to ion)
                double Fii =
                    1.0 * omega_c * nu_ii[std::make_pair(s, s2)][p] * mass_s *
                    values[ni_idx][p] *
                    (values[vi_idx2][p] / (mass_s2 * values[ni_idx2][p]) -
                     values[vi_idx][p] / (mass_s * values[ni_idx][p]));
                frictions[vi_idx][p] += Fii;
                frictions[vi_idx2][p] -= Fii;

                // Frictional heating (ion2 to ion)
                double Wii =
                    mu * Fii *
                    (values[vi_idx2][p] / (mass_s2 * values[ni_idx2][p]) -
                     values[vi_idx][p] / (mass_s * values[ni_idx][p]));

                // convert to Tnorm eV
                Wii *= mesh_length * mesh_length / (Tnorm * constants::qeomp);
                frictions[ei_idx][p] += mass_s2 * Wii;
                frictions[ei_idx2][p] += mass_s * Wii;

                // Heat exchange (ion2 to ion)
                double Qii = 3 * mu * omega_c *
                             nu_ii[std::make_pair(s, s2)][p] * mass_s *
                             values[ni_idx][p] *
                             (values[ei_idx2][p] - values[ei_idx][p]);

                frictions[ei_idx][p] += Qii;
                frictions[ei_idx2][p] -= Qii;
            }
        }
    }
}
} // namespace PENKNIFE
