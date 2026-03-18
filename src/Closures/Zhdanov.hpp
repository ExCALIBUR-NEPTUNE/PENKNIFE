#ifndef ZHDANOV_HPP
#define ZHDANOV_HPP

#include "Closure.hpp"

namespace SD = Nektar::SpatialDomains;

namespace PENKNIFE
{
// Forward declarations
class PlasmaSystem;
/**
 *
 */
class Zhdanov : public Closure
{
public:
    friend class MemoryManager<Zhdanov>;

    /// Creates an instance of this class
    static ClosureSharedPtr create(const std::weak_ptr<PlasmaSystem> &pSystem,
                                   const int spaceDim)
    {
        ClosureSharedPtr p =
            MemoryManager<Zhdanov>::AllocateSharedPtr(pSystem, spaceDim);
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

    void FillM3(double **m3)
    {
        for (int a = 0; a < Nchem; ++a)
        {
            double m_a = masses[a];
            for (int b = 0; b < Nchem; ++b)
            {
                double m_b       = masses[b];
                m3[b][a]         = G6(m_a, m_b, lambda) / pressure[b];
                m3[b][a + Nchem] = 7. * mu(m_a, m_b) * G10(m_a, m_b, lambda) /
                                   (m_a * pressure[b]);
                m3[b + Nchem][a] =
                    mu(m_a, m_b) * G10(m_a, m_b, lambda) / (T * pressure[b]);
                m3[b + Nchem][a + Nchem] =
                    m_b * G12(m_a, m_b, lambda) / (T * pressure[b]);

                m3[a][a] += G5(m_a, m_b, lambda) / pressure[a];
                m3[a][a + Nchem] += 7. * mu(m_a, m_b) * G9(m_a, m_b, lambda) /
                                    (m_a * pressure[a]);
                m3[a + Nchem][a] += 7. * mu(m_a, m_b) * G9(m_a, m_b, lambda) /
                                    (T * pressure[a]);
                m3[a + Nchem][a + Nchem] +=
                    m_a * G11(m_a, m_b, lambda) / (T * pressure[1]);
            }
        }
    }
    void FillB(double *B)
    {
        for (int a = 0; a < Nchem; ++a)
        {
            double m_a = masses[a];

            double w = 0;
            for (int b = 0; b < Nchem; ++b)
            {
                double m_b = masses[b];

                double mu_ = mu(m_a, m_b);
                w += 2.5 * mu_ * G2(m_a, m_b, lambda) * (B[b] - B[a]) / m_a;
                w += 17.5 * mu_ * mu_ * G8(m_a, m_b, lambda) *
                     (B[b + Nchem] - B[a + Nchem]) / (m_a * m_a);
            }

            B[a + Nchem] = w;
            B[a]         = 2.5 * n_bar[a] * B[a];
        }
    }

    inline void LUPDecompose(double **A, int N, double Tol, int *P)
    {

        int i, j, k, imax;
        double maxA, *ptr, absA;

        for (i = 0; i <= N; i++)
            P[i] = i; // Unit permutation matrix, P[N] initialized with N

        for (i = 0; i < N; i++)
        {
            maxA = 0.0;
            imax = i;

            for (k = i; k < N; k++)
                if ((absA = fabs(A[k][i])) > maxA)
                {
                    maxA = absA;
                    imax = k;
                }

            if (maxA < Tol)
                return; // failure, matrix is degenerate

            if (imax != i)
            {
                // pivoting P
                j       = P[i];
                P[i]    = P[imax];
                P[imax] = j;

                // pivoting rows of A
                ptr     = A[i];
                A[i]    = A[imax];
                A[imax] = ptr;

                // counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
        }
    }

    /* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
     * OUTPUT: x - solution vector of A*x=b
     */
    inline void LUPSolve(double **A, int *P, double *b, int N, double *x)
    {
        for (int i = 0; i < N; i++)
        {
            x[i] = b[P[i]];

            for (int k = 0; k < i; k++)
                x[i] -= A[i][k] * x[k];
        }

        for (int i = N - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < N; k++)
                x[i] -= A[i][k] * x[k];

            x[i] /= A[i][i];
        }
    }

    inline void Solve_qBar_rBar(double *w_bar_gradTbar, double *q_bar_r_bar)
    {
        FillB(w_bar_gradTbar);

        double **M3 = (double **)std::malloc(2 * Nchem * sizeof(double *));
        for (int i = 0; i < 2 * Nchem; i++)
            M3[i] = (double *)std::malloc(2 * Nchem * sizeof(double));

        FillM3(M3);

        int *P = (int *)std::malloc(sizeof(int) * (Nchem + 1));
        LUPDecompose(M3, Nchem, 1e-10, P);
        LUPSolve(M3, P, w_bar_gradTbar, Nchem, q_bar_r_bar);

        std::free(P);
        std::free(M3);
        for (int i = 0; i < 2 * Nchem; i++)
            std::free(M3[i]);
    }

    inline double kappa(double m_a, double m_b)
    {
        return m_a * m_b / ((m_a + m_b) * (m_a + m_b));
    }

    inline double mu(double m_a, double m_b)
    {
        return m_a * m_b / (m_a + m_b);
    }

    inline double G1(double lambda)
    {
        return -lambda;
    }
    inline double G2(double m_a, double m_b, double lambda)
    {
        return 0.6 * lambda;
    }
    inline double G3(double m_a, double m_b, double lambda)
    {
        return -2. * (1. + 0.6 * m_b / m_a) * lambda;
    }
    inline double G4(double m_a, double m_b, double lambda)
    {
        return 0.8 * lambda;
    }
    inline double G5(double m_a, double m_b, double lambda)
    {
        return -(1.3 * m_b / m_a + 1.6 + 3. * m_a / m_b) * kappa(m_a, m_b) *
               lambda;
    }
    inline double G6(double m_a, double m_b, double lambda)
    {
        return 2.7 * kappa(m_a, m_b) * lambda;
    }
    inline double G7(double m_a, double m_b, double lambda)
    {
    }
    inline double G8(double m_a, double m_b, double lambda)
    {
        return -(3. / 14.) * lambda;
    }
    inline double G9(double m_a, double m_b, double lambda)
    {
        return 0.6 * ((23. / 28.) * m_b / m_a + 8. / 7. + 3. * m_a / m_b) *
               kappa(m_a, m_b) * lambda;
    }
    inline double G10(double m_a, double m_b, double lambda)
    {
        return -(45. / 28.) * kappa(m_a, m_b) * lambda;
    }
    inline double G11(double m_a, double m_b, double lambda)
    {
        return ((433. / 280.) * (m_b * m_b) / (m_a * m_a) +
                (136. / 35.) * m_b / m_a + (459. / 35.) + 6.4 * m_a / m_b +
                3. * (m_a * m_a) / (m_b * m_b));
    }
    inline double G12(double m_a, double m_b, double lambda)
    {
        return 9.375 * kappa(m_a, m_b) * kappa(m_a, m_b) * lambda;
    }
    inline double G13(double m_a, double m_b, double lambda)
    {
        return ((18. / 35.) * m_b / m_a + 1.2) * lambda;
    }
    inline double G14(double m_a, double m_b, double lambda)
    {
        return -(24. / 35.) * lambda;
    }
    inline double G15(double m_a, double m_b, double lambda)
    {
        return -((51. / 35.) * (m_b * m_b) / (m_a * m_a) +
                 (37. / 7.) * m_b / m_a + 4.4 + 4. * m_a / m_b) *
               kappa(m_a, m_b) * lambda;
    }
    inline double G16(double m_a, double m_b, double lambda)
    {
        return (24. / 7.) * (m_b / m_a) * kappa(m_a, m_b) * lambda;
    }

    inline double S2(double m_a, double lambda)
    {
        double s_2 = 0.0;
        for (int s = 0; s < Nspec; ++s)
        {
            double m_b = masses[s];
            s_2 += mu(m_a, m_b) * G2(m_a, m_b, lambda) / m_a;
        }
        return 2.5 * s_2;
    }

    inline double S5(double m_a, double lambda)
    {
        double s_5 = 0.0;
        for (int s = 0; s < Nspec; ++s)
        {
            double m_b = masses[s];
            s_5 += G5(m_a, m_b, lambda);
        }
        return s_5;
    }

    inline double S8(double m_a, double lambda)
    {
        double s_8 = 0.0;
        for (int s = 0; s < Nspec; ++s)
        {
            double m_b = masses[s];
            s_8 += mu(m_a, m_b) * mu(m_a, m_b) * G8(m_a, m_b, lambda) /
                   (m_a * m_a);
        }
        return 17.5 * s_8;
    }

    inline double S9(double m_a, double lambda)
    {
        double s_9 = 0.0;
        for (int s = 0; s < Nspec; ++s)
        {
            double m_b = masses[s];
            s_9 += mu(m_a, m_b) * G9(m_a, m_b, lambda) / m_a;
        }
        return s_9;
    }
    inline double S11(double m_a, double lambda)
    {
        double s_11 = 0.0;
        for (int s = 0; s < Nspec; ++s)
        {
            double m_b = masses[s];
            s_11 += G11(m_a, m_b, lambda);
        }
        return s_11;
    }

    inline double c5(double nu_aa, double nu_a, double s_11, double D)
    {
        return 2.5 * nu_aa * s_11 / (nu_a * D);
    }
    inline double c6(double s_2, double s_11, double s_8, double s_9, double D)
    {
        return (s_2 * s_11 - s_8 * s_9) / D;
    }
    inline double D(double s_5, double s_11, double s_9)
    {
        s_5 *s_11 - 7. * s_9;
    }

    Zhdanov(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim);

    ~Zhdanov() = default;

    int Nchem;
    int Nspec;
    int *specs;

    std::map<std::pair<int, int>, Array<OneD, NekDouble>> nu_ii;
    std::map<int, Array<OneD, NekDouble>> nu_ei;
    Array<OneD, NekDouble> nu_ee;

    std::map<int, Array<OneD, NekDouble>> nu_i;
    Array<OneD, NekDouble> nu_e;

    NekDouble k_ci;
    NekDouble k_ce;
};

} // namespace PENKNIFE

#endif