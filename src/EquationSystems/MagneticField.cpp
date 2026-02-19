#include "MagneticField.hpp"
#include "PlasmaSystem.hpp"
#include <FieldUtils/Interpolator.h>

namespace PENKNIFE
{

MagneticField::MagneticField(const LU::SessionReaderSharedPtr &session,
                             const std::weak_ptr<PlasmaSystem> &eqn_sys,
                             Array<OneD, MR::DisContFieldSharedPtr> &pB,
                             const int dim)
    : B(pB), ndim(dim), ptsIO(session->GetComm()), m_sys(eqn_sys)
{
    npoints      = eqn_sys.lock()->GetNpoints();
    this->mag_B  = Array<OneD, NekDouble>(npoints, 0.0);
    this->b_unit = Array<OneD, Array<OneD, NekDouble>>(3);
    this->B_in   = Array<OneD, Array<OneD, NekDouble>>(3);

    for (int d = 0; d < 3; ++d)
    {
        this->b_unit[d] = Array<OneD, NekDouble>(npoints, 0.0);
        B_in[d]         = Array<OneD, NekDouble>(npoints, 0.0);
    }

    if (session->DefinesFunction("MagneticMeanField"))
    {
        if (dim == 2)
        {
            this->type = field_type::func;
            m_function = eqn_sys.lock()->GetFunction("MagneticMeanField");
        }
        else if (dim == 3)
        {
            LU::FunctionType vType =
                session->GetFunctionType("MagneticMeanField", "Bx");
            if (vType == LU::eFunctionTypeFile)
            {
                this->type = field_type::mean_file;
                filename =
                    session->GetFunctionFilename("MagneticMeanField", "Bx");
                LU::PtsFieldSharedPtr inPts2D;
                ptsIO.Import(filename, inPts2D);
                // B2D is [x,y,Bx,By,Bz]
                inPts2D->GetPts(B2D);
                unsigned int npoints_2d = inPts2D->GetNpoints();

                Array<OneD, Array<OneD, NekDouble>> pts(6);
                for (int i = 0; i < 3; ++i)
                {
                    pts[i]     = Array<OneD, NekDouble>(npoints);
                    pts[i + 3] = B_in[i];
                }
                B[0]->GetCoords(pts[0], pts[1], pts[2]);

                this->B3D = Array<OneD, Array<OneD, NekDouble>>(6);
                unsigned int increments = 64;

                for (int d = 0; d < 6; ++d)
                {
                    B3D[d] =
                        Array<OneD, NekDouble>(increments * npoints_2d, 0.0);
                }

                for (int i = 0; i < increments; ++i)
                {
                    NekDouble theta = i * 2 * M_PI / increments;
                    for (int j = 0; j < npoints_2d; ++j)
                    {
                        int p     = i * npoints_2d + j;
                        B3D[0][p] = cos(theta) * B2D[0][j];
                        B3D[1][p] = B2D[1][j];
                        B3D[2][p] = sin(theta) * B2D[0][j];
                    }
                }
                inPts  = MemoryManager<LU::PtsField>::AllocateSharedPtr(3, B3D);
                outPts = MemoryManager<LU::PtsField>::AllocateSharedPtr(
                    3, Bstring, pts);
                interp =
                    FieldUtils::Interpolator<std::vector<MR::ExpListSharedPtr>>(
                        LU::eShepard);
                interp.CalcWeights(inPts, outPts);
            }
            else if (vType == LU::eFunctionTypeExpression)
            {
                this->type = field_type::mean_exp;

                Bxfunc = session->GetFunction("MagneticMeanField", "Bx");
                Byfunc = session->GetFunction("MagneticMeanField", "By");
                Bzfunc = session->GetFunction("MagneticMeanField", "Bz");

                // Get the coordinates (assuming all fields have the same
                // discretisation)
                Array<OneD, NekDouble> x0(npoints);
                x1 = Array<OneD, NekDouble>(npoints);
                Array<OneD, NekDouble> x2(npoints);
                B[0]->GetCoords(x0, x1, x2);
                r   = Array<OneD, NekDouble>(npoints);
                phi = Array<OneD, NekDouble>(npoints);

                Br   = Array<OneD, NekDouble>(npoints);
                Bphi = Array<OneD, NekDouble>(npoints);

                for (int q = 0; q < npoints; ++q)
                {
                    r[q]   = std::sqrt(x0[q] * x0[q] + x2[q] * x2[q]);
                    phi[q] = std::atan2(x2[q], x0[q]);
                }
            }
        }
    }
    else if (session->DefinesFunction("MagneticField"))
    {
        this->type = field_type::func;
        m_function = eqn_sys.lock()->GetFunction("MagneticField");
    }
    else if (session->DefinesFunction("CurrentDensity"))
    {
        this->type = field_type::solve;
        m_function = eqn_sys.lock()->GetFunction("CurrentDensity");
    }
}
/**
 * @brief Read the magnetic field from file.
 */
void MagneticField::Read(NekDouble time)
{
    int d;

    if (this->type == field_type::mean_file)
    {
        LU::PtsFieldSharedPtr inPts2D;
        ptsIO.Import(filename, inPts2D);
        // B2D is [x,y,Bx,By,Bz]
        inPts2D->GetPts(B2D);
        unsigned int npoints_2d = inPts2D->GetNpoints();
        // B3D is [x,y,z,Bx,By,Bz]
        for (int i = 0; i < increments; ++i)
        {
            NekDouble theta = i * 2 * M_PI / increments;
            for (int j = 0; j < npoints_2d; ++j)
            {
                int p     = i * npoints_2d + j;
                B3D[3][p] = cos(theta) * B2D[2][j] - sin(theta) * B2D[4][j];
                B3D[4][p] = B2D[3][j];
                B3D[5][p] = cos(theta) * B2D[4][j] + sin(theta) * B2D[2][j];
            }
        }

        interp.Interpolate(inPts, outPts);
    }
    else if (this->type == field_type::mean_exp)
    {
        Bxfunc->Evaluate(r, x1, phi, time, Br);
        Byfunc->Evaluate(r, x1, phi, time, B_in[1]);
        Bzfunc->Evaluate(r, x1, phi, time, Bphi);

        for (int q = 0; q < npoints; ++q)
        {
            B_in[0][q] = cos(phi[q]) * Br[q] - sin(phi[q]) * Bphi[q];
            B_in[2][q] = cos(phi[q]) * Bphi[q] + sin(phi[q]) * Br[q];
        }
    }

    else if (this->type == field_type::func)
    {
        m_function->Evaluate(Bstring, B_in, time);
    }
    else
    {
        ASSERTL0(false, "Unknown B field type");
    }

    for (d = 0; d < 3; ++d)
    {
        this->B[d]->UpdatePhys() = B_in[d];
    }
}

void MagneticField::Solve(Array<OneD, Array<OneD, NekDouble>> &J)
{
    m_function->Evaluate("Bz", B_in[2]);
    Array<OneD, NekDouble> A[3];
    Array<OneD, NekDouble> Acoeff[3];
    int ncoeffs = m_sys.lock()->GetNcoeffs();

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    factors[StdRegions::eFactorTau]    = 1.0;

    for (int d = 0; d < 3; ++d)
    {
        Acoeff[d] = Array<OneD, NekDouble>(ncoeffs, 0.0);
        A[d]      = Array<OneD, NekDouble>(npoints, 0.0);
        B[d]->HelmSolve(B_in[d], Acoeff[d], factors);
        B[d]->BwdTrans(Acoeff[d], A[d]);
    }

    if (ndim == 3)
    {
        Array<OneD, NekDouble> tmp1(npoints, 0.0);
        Array<OneD, NekDouble> tmp2(npoints, 0.0);
        B[0]->PhysDeriv(1, A[2], tmp1);
        B[0]->PhysDeriv(2, A[1], tmp2);
        Vmath::Vsub(npoints, tmp1, 1, tmp2, 1, B[0]->UpdatePhys(), 1);
        B[1]->PhysDeriv(2, A[0], tmp1);
        B[1]->PhysDeriv(0, A[2], tmp2);
        Vmath::Vsub(npoints, tmp1, 1, tmp2, 1, B[1]->UpdatePhys(), 1);
        B[2]->PhysDeriv(0, A[1], tmp1);
        B[2]->PhysDeriv(1, A[0], tmp2);
        Vmath::Vsub(npoints, tmp1, 1, tmp2, 1, B[2]->UpdatePhys(), 1);
    }
    else if (ndim == 2)
    {
        B[0]->PhysDeriv(1, A[2], B[0]->UpdatePhys());
        B[1]->PhysDeriv(0, A[2], B[1]->UpdatePhys());
        Vmath::Neg(npoints, B[1]->UpdatePhys(), 1);
    }
}

void MagneticField::Update(NekDouble time)
{
    if (this->type == field_type::solve)
    {
        Solve();
    }
    else
    {
        Read(time);
    }
    for (int d = 0; d < 3; ++d)
    {
        this->B[d]->FwdTrans(this->B[d]->GetPhys(), this->B[d]->UpdateCoeffs());
        Vmath::Vvtvp(npoints, this->B[d]->GetPhys(), 1, this->B[d]->GetPhys(),
                     1, this->mag_B, 1, this->mag_B, 1);
    }

    for (int d = 0; d < 3; d++)
    {
        for (int k = 0; k < npoints; ++k)
        {
            this->b_unit[d][k] =
                (this->mag_B[k] > 0)
                    ? this->B[d]->GetPhys()[k] / std::sqrt(this->mag_B[k])
                    : 0.0;
        }
    }
}

} // namespace PENKNIFE