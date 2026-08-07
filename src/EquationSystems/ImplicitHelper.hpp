///////////////////////////////////////////////////////////////////////////////
//
// File DriftWaveSystem.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Unsteady driftwave solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_IMPLICITHELPER_HPP
#define NEKTAR_IMPLICITHELPER_HPP

#include <LibUtilities/LinearAlgebra/NekNonlinSysIter.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

/**
 * @brief An equation system for the drift-wave solver.
 */
class ImplicitHelper
{
public:
    ImplicitHelper(LibUtilities::SessionReaderSharedPtr session,
                   Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                   LibUtilities::TimeIntegrationSchemeOperators &ode,
                   int nFields)
        : m_session(session), m_fields(fields), m_ode(ode), m_nFields(nFields)
    {
        m_comm = session->GetComm()->GetSpaceComm();
    }

    /// Default destructor.
    ~ImplicitHelper() = default;

public:
    void ImplicitTimeInt(
        const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
        Array<OneD, Array<OneD, NekDouble>> &outpnt, const NekDouble time,
        const NekDouble lambda)
    {
        m_TimeIntegLambda    = lambda;
        m_bndEvaluateTime    = time;
        unsigned int npoints = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> inarray(m_nFields * npoints);
        Array<OneD, NekDouble> outarray(m_nFields * npoints);
        Array<OneD, NekDouble> tmp;

        for (int i = 0; i < m_nFields; ++i)
        {
            int noffset = i * npoints;
            Vmath::Vcopy(npoints, inpnts[i], 1, tmp = inarray + noffset, 1);
        }

        ImplicitTimeInt1D(inarray, outarray);

        for (int i = 0; i < m_nFields; ++i)
        {
            int noffset = i * npoints;
            Vmath::Vcopy(npoints, outarray + noffset, 1, outpnt[i], 1);
        }

    }

    void InitialiseNonlinSysSolver(bool coeff = false)
    {
        int ntotal = m_nFields * m_fields[0]->GetNpoints();

        // Create the key to hold settings for nonlin solver
        LibUtilities::NekSysKey key = LibUtilities::NekSysKey();

        // Load required LinSys parameters:
        m_session->LoadParameter("NekLinSysMaxIterations",
                                 key.m_NekLinSysMaxIterations, 30);
        m_session->LoadParameter("LinSysMaxStorage", key.m_LinSysMaxStorage,
                                 30);
        m_session->LoadParameter("LinSysRelativeTolInNonlin",
                                 key.m_NekLinSysTolerance, 5.0E-2);
        m_session->LoadParameter("GMRESMaxHessMatBand",
                                 key.m_KrylovMaxHessMatBand, 31);

        // Load required NonLinSys parameters:
        m_session->LoadParameter("JacobiFreeEps", m_jacobiFreeEps, 5.0E-8);
        m_session->LoadParameter("NekNonlinSysMaxIterations",
                                 key.m_NekNonlinSysMaxIterations, 10);
        m_session->LoadParameter("NewtonRelativeIteTol",
                                 key.m_NekNonLinSysTolerance, 1.0E-12);
        WARNINGL0(!m_session->DefinesParameter("NewtonAbsoluteIteTol"),
                  "Please specify NewtonRelativeIteTol instead of "
                  "NewtonAbsoluteIteTol in XML session file");
        m_session->LoadParameter("NonlinIterTolRelativeL2",
                                 key.m_NonlinIterTolRelativeL2, 1.0E-3);
        m_session->LoadSolverInfo("LinSysIterSolverTypeInNonlin",
                                  key.m_LinSysIterSolverTypeInNonlin, "GMRES");

        // Set up operators
        LibUtilities::NekSysOperators nekSysOp;
        if (coeff)
        {
            nekSysOp.DefineNekSysResEval(
                &ImplicitHelper::NonlinSysEvaluatorCoeff1D, this);
            nekSysOp.DefineNekSysLhsEval(
                &ImplicitHelper::MatrixMultiplyMatrixFreeCoeff, this);
        }
        else
        {
            nekSysOp.DefineNekSysResEval(&ImplicitHelper::NonlinSysEvaluator1D,
                                         this);
            nekSysOp.DefineNekSysLhsEval(
                &ImplicitHelper::MatrixMultiplyMatrixFree, this);
        }
        nekSysOp.DefineNekSysPrecon(&ImplicitHelper::DoNullPrecon, this);

        // Initialize non-linear system
        m_nonlinsol = LibUtilities::GetNekNonlinSysIterFactory().CreateInstance(
            "Newton", m_session, m_comm->GetRowComm(), ntotal, key);
        m_nonlinsol->SetSysOperators(nekSysOp);
    }

protected:
    // Implicit solver parameters
    int m_nFields               = 0;
    int m_TotNewtonIts          = 0;
    int m_TotLinIts             = 0;
    int m_TotImpStages          = 0;
    NekDouble m_jacobiFreeEps   = 5.0E-08;
    NekDouble m_bndEvaluateTime = 0.0;
    NekDouble m_TimeIntegLambda = 0.0;
    NekDouble m_inArrayNorm     = -1.0;

    LibUtilities::SessionReaderSharedPtr m_session;
    LibUtilities::NekNonlinSysIterSharedPtr m_nonlinsol;
    LibUtilities::CommSharedPtr m_comm;
    LibUtilities::TimeIntegrationSchemeOperators &m_ode;

    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

    void ImplicitTimeInt1D(const Array<OneD, const NekDouble> &inarray,
                           Array<OneD, NekDouble> &out)
    {
        CalcRefValues(inarray);
        m_nonlinsol->SetRhsMagnitude(m_inArrayNorm);
        m_TotNewtonIts +=
            m_nonlinsol->SolveSystem(inarray.size(), inarray, out, 0);
        m_TotLinIts += m_nonlinsol->GetNtotLinSysIts();
        m_TotImpStages++;
    }

    void CalcRefValues(const Array<OneD, const NekDouble> &inarray)
    {
        unsigned int npoints = m_fields[0]->GetNpoints();

        Array<OneD, NekDouble> magnitdEstimat(m_nFields, 0.0);

        for (int i = 0; i < m_nFields; ++i)
        {
            int offset = i * npoints;
            magnitdEstimat[i] =
                Vmath::Dot(npoints, inarray + offset, inarray + offset);
        }
        m_comm->GetSpaceComm()->AllReduce(magnitdEstimat,
                                          Nektar::LibUtilities::ReduceSum);

        m_inArrayNorm = 0.0;
        for (int i = 0; i < m_nFields; ++i)
        {
            m_inArrayNorm += magnitdEstimat[i];
        }
    }

    void NonlinSysEvaluator1D(const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &out,
                              [[maybe_unused]] const bool &flag)
    {
        unsigned int npoints = m_fields[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble>> in2D(m_nFields);
        Array<OneD, Array<OneD, NekDouble>> out2D(m_nFields);
        for (int i = 0; i < m_nFields; ++i)
        {
            int offset = i * npoints;
            in2D[i]    = inarray + offset;
            out2D[i]   = out + offset;
        }
        NonlinSysEvaluator(in2D, out2D);
    }

    void NonlinSysEvaluator(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &out)
    {
        unsigned int npoints = m_fields[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble>> inpnts(m_nFields);
        for (int i = 0; i < m_nFields; ++i)
        {
            inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        m_ode.DoProjection(inarray, inpnts, m_bndEvaluateTime);
        m_ode.DoOdeRhs(inpnts, out, m_bndEvaluateTime);

        for (int i = 0; i < m_nFields; ++i)
        {
            Vmath::Svtvp(npoints, -m_TimeIntegLambda, out[i], 1, inarray[i], 1,
                         out[i], 1);
            Vmath::Vsub(npoints, out[i], 1,
                        m_nonlinsol->GetRefSourceVec() + i * npoints, 1, out[i],
                        1);
        }
    }
    void NonlinSysEvaluatorCoeff1D(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &out,
                                   [[maybe_unused]] const bool &flag)
    {
        unsigned int ncoeffs = m_fields[0]->GetNcoeffs();
        Array<OneD, Array<OneD, NekDouble>> in2D(m_nFields);
        Array<OneD, Array<OneD, NekDouble>> out2D(m_nFields);
        for (int i = 0; i < m_nFields; ++i)
        {
            int offset = i * ncoeffs;
            in2D[i]    = inarray + offset;
            out2D[i]   = out + offset;
        }

        NonlinSysEvaluatorCoeff(in2D, out2D);
    }

    void NonlinSysEvaluatorCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &out)
    {
        unsigned int ncoeffs = m_fields[0]->GetNcoeffs();
        unsigned int npoints = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble>> inpnts(m_nFields);

        for (int i = 0; i < m_nFields; ++i)
        {
            inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }

        m_ode.DoProjection(inpnts, inpnts, m_bndEvaluateTime);
        m_ode.DoOdeRhs(inpnts, out, m_bndEvaluateTime);
        for (int i = 0; i < m_nFields; ++i)
        {
            Vmath::Svtvp(ncoeffs, -m_TimeIntegLambda, out[i], 1, inarray[i], 1,
                         out[i], 1);
            Vmath::Vsub(ncoeffs, out[i], 1,
                        m_nonlinsol->GetRefSourceVec() + i * ncoeffs, 1, out[i],
                        1);
        }
    }

    void MatrixMultiplyMatrixFree(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out,
                                  [[maybe_unused]] const bool &flag)
    {
        const Array<OneD, const NekDouble> solref =
            m_nonlinsol->GetRefSolution();
        const Array<OneD, const NekDouble> resref =
            m_nonlinsol->GetRefResidual();

        unsigned int ntotal   = inarray.size();
        NekDouble magninarray = Vmath::Dot(ntotal, inarray, inarray);
        m_comm->GetSpaceComm()->AllReduce(magninarray,
                                          Nektar::LibUtilities::ReduceSum);
        NekDouble eps =
            m_jacobiFreeEps * sqrt((sqrt(m_inArrayNorm) + 1.0) / magninarray);

        Array<OneD, NekDouble> solplus{ntotal};
        Array<OneD, NekDouble> resplus{ntotal};

        Vmath::Svtvp(ntotal, eps, inarray, 1, solref, 1, solplus, 1);
        NonlinSysEvaluator1D(solplus, resplus, flag);
        Vmath::Vsub(ntotal, resplus, 1, resref, 1, out, 1);
        Vmath::Smul(ntotal, 1.0 / eps, out, 1, out, 1);
    }

    void MatrixMultiplyMatrixFreeCoeff(
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble> &out, [[maybe_unused]] const bool &flag)
    {
        const Array<OneD, const NekDouble> solref =
            m_nonlinsol->GetRefSolution();
        const Array<OneD, const NekDouble> resref =
            m_nonlinsol->GetRefResidual();

        unsigned int ntotal   = inarray.size();
        NekDouble magninarray = Vmath::Dot(ntotal, inarray, inarray);
        m_comm->GetSpaceComm()->AllReduce(magninarray,
                                          Nektar::LibUtilities::ReduceSum);
        NekDouble eps =
            m_jacobiFreeEps * sqrt((sqrt(m_inArrayNorm) + 1.0) / magninarray);

        Array<OneD, NekDouble> solplus{ntotal};
        Array<OneD, NekDouble> resplus{ntotal};

        Vmath::Svtvp(ntotal, eps, inarray, 1, solref, 1, solplus, 1);
        NonlinSysEvaluatorCoeff1D(solplus, resplus, flag);
        Vmath::Vsub(ntotal, resplus, 1, resref, 1, out, 1);
        Vmath::Smul(ntotal, 1.0 / eps, out, 1, out, 1);
    }

    void DoNullPrecon(const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray, const bool &flag)
    {
        Vmath::Vcopy(inarray.size(), inarray, 1, outarray, 1);
    }
};

} // namespace Nektar

#endif