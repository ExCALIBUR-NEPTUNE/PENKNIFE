#include "Zhdanov.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

namespace PENKNIFE
{
std::string Zhdanov::className = GetClosureFactory().RegisterCreatorFunction(
    "Zhdanov", Zhdanov::create, "Zhdanov closure system");

Zhdanov::Zhdanov(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim)
    : Closure(pSystem, spaceDim)
{
    Nchem        = m_system.lock()->GetChem().size();
    this->specs  = (int *)std::malloc(sizeof(int) * (Nchem + 1));
    this->mass   = (double *)std::malloc(sizeof(double) * (Nchem + 1));
    this->charge = (double **)std::malloc(sizeof(double *) * Nchem);
    Nspec        = m_system.lock()->GetIons().size();
    this->n_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));
    this->v_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));
    this->e_idx  = (int *)std::malloc(sizeof(int) * (Nspec + 1));

    int cnt = 0;
    for (const auto &[c, cv] : m_system.lock()->GetChem())
    {
        specs[c]  = 0;
        charge[c] = (double *)std::malloc(sizeof(double) * cv.specs.size());
        mass[c]   = cv.mass;
        for (int s = 0; s < cv.specs.size(); ++s)
        {
            specs[c]++;
            charge[c][s] = cv.specs[s].charge;
            n_idx[cnt]   = cv.specs[s].fields.at(field_to_index.at("n"));
            v_idx[cnt]   = cv.specs[s].fields.at(field_to_index.at("v"));
            e_idx[cnt]   = cv.specs[s].fields.at(field_to_index.at("e"));
            idx[c][s]    = cnt++;
        }
    }

    // Electrons included as "element" with single charge state
    this->specs[Nchem]     = 1;
    this->mass[Nchem]      = constants::m_e_m_p;
    this->charge[Nchem][0] = -1;
    Nchem++;
    e_idx[Nspec]     = ee_idx;
    n_idx[Nchem - 1] = m_system.lock()->n_indep_fields;
    v_idx[Nchem - 1] = m_system.lock()->n_indep_fields + 1;
    Nspec++;

    // Set up maps from indices to collision parameters
    for (int a = 0; a < Nchem; ++a)
    {
        this->nu_a[a] = Array<OneD, NekDouble>(this->n_pts, 0.0);

        for (int b = 0; b < Nchem; ++b)
        {
            if (b > a)
                break;
            this->lambda_ab[{a, b}] = Array<OneD, NekDouble>(this->n_pts, 0.0);

            for (int z = 0; z < specs[a]; ++z)
            {
                for (int y = 0; y < specs[b]; ++y)
                {
                    if (b == a && y > z)
                        break;
                    int idx1 = idx[a][z];
                    int idx2 = idx[b][y];
                    this->lambda_aZbY[{idx1, idx2}] =
                        Array<OneD, NekDouble>(this->n_pts);
                }
            }
        }
    }
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

inline double CoulombLog(double Nnorm, double Tnorm, double ni1, double ni2,
                         double Ti1, double Ti2, double A1, double A2,
                         double Z1, double Z2)
{
}

void Zhdanov::CalcLambdas(const Array<OneD, Array<OneD, NekDouble>> &in_arr,
                          const Array<OneD, NekDouble> &ne)
{
    for (int p = 0; p < this->n_pts; ++p)
    {
        for (int a = 0; a < Nchem; ++a)
        {
            double A   = mass[a];
            nu_a[a][p] = 0.0;
            for (int b = 0; b < Nchem; ++b)
            {
                if (b > a)
                    break;
                double B   = mass[b];
                double mu_ = mu(A, B);
                for (int z = 0; z < specs[a]; ++z)
                {
                    double Z = charge[a][z];
                    for (int y = 0; y < specs[b]; ++y)
                    {
                        if (b == a && y > z)
                            break;
                        int idx1 = idx[a][z];
                        int idx2 = idx[a][z];

                        double Y           = charge[b][y];
                        double coulomb_log = CoulombLog(
                            in_arr[n_idx[idx1]][p], in_arr[n_idx[idx2]][p],
                            in_arr[e_idx[idx1]][p], in_arr[e_idx[idx2]][p], A,
                            B, Z, Y);

                        const double v1sq =
                            2 * Tnorm * in_arr[e_idx[idx1]][p] / A;
                        const double v2sq =
                            2 * Tnorm * in_arr[e_idx[idx2]][p] / B;

                        double lambda =
                            Z * Z * Y * Y * in_arr[n_idx[idx1]][p] *
                            in_arr[n_idx[idx2]][p] * coulomb_log /
                            (3 * mu_ * pow(M_PI * (v1sq + v2sq), 1.5) *
                             pow(constants::epsilon_0, 2));

                        lambda *= (sqrt(constants::m_p) / constants::c) *
                                  Nnorm * Nnorm / 1e12;
                        lambda *= (pow(mesh_length, 3) / omega_c);

                        lambda_aZbY[{idx1, idx2}][p] = lambda;
                        lambda_ab[{a, b}][p] += lambda;
                        if (b == a && y < z)
                            lambda_ab[{a, b}][p] += lambda;
                    }
                }
                nu_a[a][p] += lambda_ab[{a, b}][p];
                if (b < a)
                    nu_a[b][p] += lambda_ab[{a, b}][p];
            }
            double n_tot = 0.0;
            for (int z = 0; z < specs[a]; ++z)
            {
                int idx = this->idx[a][z];
                n_tot += in_arr[n_idx[idx]][p];
            }
            nu_a[a][p] /= (n_tot * mass[a]);
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
    // Append electron data to ion data
    Array<OneD, Array<OneD, NekDouble>> values(vals.size() + 2);
    for (int v = 0; v < vals.size(); v++)
        values[v] = vals[v];
    values[vals.size()]     = ne;
    values[vals.size() + 1] = ve;
    CalcLambdas(values, ne);

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

        Solve_qBar_rBar(p, w_bar_gradTbar, q_bar_r_bar, n_bar, p_bar, T);
        double *q_bar = q_bar_r_bar;
        double *r_bar = q_bar_r_bar + Nchem;

        double *r = (double *)std::malloc(sizeof(double) * Nspec);
        double *q = (double *)std::malloc(sizeof(double) * Nspec);

        for (int c = 0; c < Nchem; ++c)
        {
            double s2, s5, s8, s9, s11;
            Calc_S_coeffs(p, c, &s2, &s5, &s8, &s9, &s11);

            double nu_aa = 2 * lambda_ab[{c, c}][p] / (n_bar[c] * mass[c]);
            double nu_a  = this->nu_a[c][p];
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
                q[i] *= pressure;
                r[i] *= pressure;
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

                // Frictional heat exchange
                Q += R * (mu_ / mass[a]) *
                     (values[v_idx[b]][p] / (values[n_idx[b]][p] * mass[b]) -
                      values[v_idx[a]][p] / (values[n_idx[a]][p] * mass[a]));
                // Thermal heat exchange
                Q += 3 * lambda * (values[e_idx[b]][p] - values[e_idx[a]][p]) /
                     (mass[a] + mass[b]);
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
