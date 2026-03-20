#include "Zhdanov.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

namespace PENKNIFE
{
std::string Zhdanov::className = GetClosureFactory().RegisterCreatorFunction(
    "Zhdanov", Zhdanov::create, "Zhdanov closure system");

Zhdanov::Zhdanov(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim)
    : Closure(pSystem, spaceDim)
{

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
    Nchem        = m_system.lock()->GetChem().size();
    this->specs  = (int *)std::malloc(sizeof(int) * (Nchem + 1));
    this->mass   = (double *)std::malloc(sizeof(double) * (Nchem + 1));
    this->charge = (double **)std::malloc(sizeof(double *) * Nchem);
    Nspec        = m_system.lock()->GetIons().size();
    this->n_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));
    this->v_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));
    this->e_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));

    int speccnt  = 0;
    int fieldcnt = 0;
    for (const auto &[c, cv] : m_system.lock()->GetChem())
    {
        specs[c]  = 0;
        charge[c] = (double *)std::malloc(sizeof(double) * cv.specs.size());
        mass[c]   = cv.mass;
        for (int s = 0; s < cv.specs.size(); ++s)
        {
            specs[c]++;
            charge[c][s]   = cv.specs[s].charge;
            n_idx[speccnt] = cv.specs[s].fields.at(field_to_index.at("n"));
            fieldcnt++;
            v_idx[speccnt] = cv.specs[s].fields.at(field_to_index.at("v"));
            fieldcnt++;
            e_idx[speccnt] = cv.specs[s].fields.at(field_to_index.at("e"));
            fieldcnt++;
            speccnt++;
        }
    }

    this->specs[Nchem]     = 1;
    this->mass[Nchem]      = constants::m_e_m_p;
    this->charge[Nchem][0] = -1;
    Nchem++;
    e_idx[Nspec]     = ee_idx;
    n_idx[Nchem - 1] = m_system.lock()->n_indep_fields;
    v_idx[Nchem - 1] = m_system.lock()->n_indep_fields + 1;
    Nspec++;
}

Zhdanov::~Zhdanov()
{
    std::free(specs);
    std::free(mass);
    for (int c = 0; c < Nchem; ++c)
    {
        std::free(charge[c]);
    }
    std::free(charge);
    std::free(n_idx);
    std::free(v_idx);
    std::free(e_idx);
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

void Zhdanov::CalcCollisionFrequencies(
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

void Zhdanov::CalcLambdas(const Array<OneD, Array<OneD, NekDouble>> &in_arr,
                          const Array<OneD, NekDouble> &ne)
{
    for (int p = 0; p < this->n_pts; ++p)
    {
        for (int a = 0; a < Nchem; ++a)
        {
            double A = mass[a];
            for (int b = 0; b < Nchem; ++b)
            {
                if (b > a)
                    break;
                double B          = mass[b];
                double mu_        = mu(A, B);
                double lambda_tot = 0.0;
                double n_tot      = 0.0;
                for (int z = 0; z < specs[a]; ++z)
                {
                    double Z = charge[a][z];
                    for (int y = 0; y < specs[b]; ++y)
                    {
                        if (b == a && y > z)
                            break;
                        double Y           = charge[b][y];
                        double coulomb_log = CoulombLog_ii(
                            Nnorm, in_arr[n_idx[z]][p], in_arr[n_idx[y]][p],
                            in_arr[e_idx[z]][p], in_arr[e_idx[y]][p], A, B, Z,
                            Y);

                        const double v1sq = 2 * in_arr[e_idx[z]][p] *
                                            constants::e /
                                            (A * constants::m_p_si);
                        const double v2sq = 2 * in_arr[e_idx[y]][p] *
                                            constants::e /
                                            (B * constants::m_p_si);

                        double lambda = Z * Z * Y * Y * sqrt(mu_) *
                                        pow(constants::e, 4) *
                                        in_arr[n_idx[z]][p] *
                                        in_arr[n_idx[y]][p] * coulomb_log /
                                        (3 * pow(M_PI * (v1sq + v2sq), 1.5) *
                                         pow(constants::epsilon_0_si, 2));
                        // lambda /= mesh_length;

                        lambda_aZbY[{z, y}][p] = lambda;
                        lambda_tot += in_arr[n_idx[z]][p] * lambda;
                        n_tot += in_arr[n_idx[z]][p];
                    }
                }
                lambda_ab[{a, b}][p] = lambda_tot / n_tot;
            }
        }
    }
}

void Zhdanov::v_EvaluateClosure(
    const Array<OneD, Array<OneD, NekDouble>> &vals,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
    Array<OneD, Array<OneD, NekDouble>> &frictions,
    const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve)
{
    CalcCollisionFrequencies(vals, ne);

    Array<OneD, Array<OneD, NekDouble>> values(vals.size() + 2);
    for (int v = 0; v < vals.size(); v++)
        values[v] = vals[v];
    values[vals.size()]     = ne;
    values[vals.size() + 1] = ve;

    for (int p = 0; p < this->n_pts; ++p)
    {
        double T         = 0;
        double total_n   = 0;
        double com_v_par = 0.0;
        double tot_mass  = 0.0;

        for (int s = 0; s < Nspec; ++s)
        {
            com_v_par += values[v_idx[s]][p];
            tot_mass += mass[s] * values[n_idx[s]][p];

            T += values[e_idx[s]][p];
            total_n += values[n_idx[s]][p];
        }

        com_v_par /= tot_mass;
        T /= total_n;
        double *w_bar_gradTbar =
            (double *)std::malloc(sizeof(double) * 2 * Nchem);

        double *w_bar    = w_bar_gradTbar;
        double *gradTbar = w_bar_gradTbar + Nchem;

        double *n_bar    = (double *)std::malloc(sizeof(double) * Nchem);
        double *p_bar    = (double *)std::malloc(sizeof(double) * Nchem);
        double *Z_sq_bar = (double *)std::malloc(sizeof(double) * Nchem);
        double *w        = (double *)std::malloc(sizeof(double) * Nspec);

        for (int c = 0; c < Nchem; ++c)
        {
            n_bar[c]     = 0.0;
            double Z_bar = 0.0;
            Z_sq_bar[c]  = 0.0;
            w_bar[c]     = 0.0;
            gradTbar[c]  = 0.0;
            int i        = 0;
            for (int s = 0; s < specs[c]; ++s)
            {
                w[i] = values[v_idx[s]][p] / (mass[c] * values[n_idx[s]][p]) -
                       com_v_par;
                double Z = charge[c][s];
                n_bar[c] += values[n_idx[s]][p];
                Z_bar += Z * values[n_idx[s]][p];
                Z_sq_bar[c] += Z * Z * values[n_idx[s]][p];
                w_bar[c] += w[i] * values[n_idx[s]][p] * Z * Z;
                for (int d = 0; d < m_spacedim; ++d)
                {
                    gradTbar[c] +=
                        Z * values[n_idx[s]][p] * grads[d][e_idx[s]][p];
                }
                p_bar[c] += values[n_idx[s]][p] * values[e_idx[s]][p];
                i++;
            }
            Z_bar /= n_bar[c];
            Z_sq_bar[c] /= n_bar[c];
            w_bar[c] /= (n_bar[c] * Z_sq_bar[c]);
            gradTbar[c] /= (n_bar[c] * Z_bar);
        }
        double *q_bar_r_bar = (double *)std::malloc(sizeof(double) * 2 * Nchem);

        Solve_qBar_rBar(w_bar_gradTbar, q_bar_r_bar, n_bar, p_bar, T);
        double *q_bar = q_bar_r_bar;
        double *r_bar = q_bar_r_bar + Nchem;

        double *r = (double *)std::malloc(sizeof(double) * Nspec);
        double *q = (double *)std::malloc(sizeof(double) * Nspec);

        for (int c = 0; c < Nchem; ++c)
        {
            double s2, s5, s8, s9, s11;
            Calc_S_coeffs(p, c, &s2, &s5, &s8, &s9, &s11);

            double nu_aa = nu_ii[{c, c}][p];
            double nu_a  = nu_i[c][p];
            double D_    = D(s5, s11, s9);
            double c_5   = c5(nu_aa, nu_a, s11, D_);
            double c_6   = c6(s2, s11, s8, s9, D_);

            double gradTcoeff = n_bar[c] * c_5 * nu_aa / nu_a;

            int i = 0;
            for (int s = 0; s < specs[c]; ++s)
            {
                double T = values[e_idx[s]][p];
                double Z = charge[c][s];
                double w =
                    values[v_idx[s]][p] / (mass[c] * values[n_idx[s]][p]) -
                    com_v_par;

                double gradTdiff = 0.0;
                for (int d = 0; d < m_spacedim; ++d)
                {
                    gradTdiff += b_unit[d][p] * grads[d][e_idx[s]][p];
                }
                gradTdiff = (Z_sq_bar[c] / Z * Z) * gradTdiff - gradTbar[c];

                q[i] = gradTcoeff * gradTdiff + c_6 * (w - w_bar[c]) +
                       q_bar[c] / p_bar[c];
                r[i] = T / (mass[c] * s11) *
                           (-7. * s9 * gradTcoeff * gradTdiff -
                            (s8 + 7. * s9 * c_6) * (w - w_bar[c])) +
                       r_bar[c] / p_bar[c];

                double pressure = (values[n_idx[s]][p] * T);
                q[i] /= pressure;
                r[i] /= pressure;
                for (int d = 0; d < m_spacedim; ++d)
                {
                    fluxes[d][e_idx[s]][p] = b_unit[d][p] * q[i];
                }
                i++;
            }
        }
        std::free(p_bar);
        std::free(n_bar);
        std::free(Z_sq_bar);
        std::free(w_bar_gradTbar);

        for (int a = 0; a < Nspec; ++a)
        {
            double R = 0;
            double Q = 0;
            for (int b = 0; b < Nspec; ++b)
            {
                double lambda = this->lambda_aZbY[{a, b}][p];
                double mu_    = mu(mass[a], mass[b]);

                double g1 = G1(lambda);
                double g2 = G2(mass[a], mass[b], lambda);
                double g8 = G8(mass[a], mass[b], lambda);
                double R1 = g1 * (w[a] - w[b]);
                double R2 = g2 * (mu_ / T) *
                            (q[a] / (values[n_idx[a]][p] * mass[a]) -
                             q[b] / (values[n_idx[b]][p] * mass[b]));
                double R3 = g8 * (mu_ * mu_ / (T * T)) *
                            (r[a] / (values[n_idx[a]][p] * mass[a]) -
                             r[b] / (values[n_idx[b]][p] * mass[b]));
                R += (R1 + R2 + R3);
                Q += R * mu_ * (values[v_idx[b]][p] - values[v_idx[a]][p]);
                Q += mu_ * nu_ii[{a, b}][p] * values[n_idx[a]][p] *
                     (values[e_idx[b]][p] - values[e_idx[a]][p]);
            }

            frictions[v_idx[a]][p] = R;
            frictions[e_idx[a]][p] = Q;
        }
        std::free(q);
        std::free(r);
        std::free(w);
    }
}
} // namespace PENKNIFE
