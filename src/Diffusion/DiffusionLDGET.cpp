#include <iomanip>
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>

#include "DiffusionLDGET.hpp"

namespace PENKNIFE
{
std::string DiffusionLDGET::type =
    GetDiffusionFactory().RegisterCreatorFunction(
        "LDGET", DiffusionLDGET::create, "Local Discontinuous Galkerin");

DiffusionLDGET::DiffusionLDGET()
{
}

void DiffusionLDGET::v_InitObject(
    LibUtilities::SessionReaderSharedPtr pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    m_session = pSession;

    m_session->LoadSolverInfo("ShockCaptureType", m_shockCaptureType, "Off");

    // Set up penalty term for LDG
    m_session->LoadParameter("LDGc11", m_C11, 1.0);

    // Setting up the normals
    std::size_t nDim      = pFields[0]->GetCoordim(0);
    std::size_t nTracePts = pFields[0]->GetTrace()->GetTotPoints();

    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>{nDim};
    for (std::size_t i = 0; i < nDim; ++i)
    {
        m_traceNormals[i] = Array<OneD, NekDouble>{nTracePts};
    }
    pFields[0]->GetTrace()->GetNormals(m_traceNormals);
}

void DiffusionLDGET::v_Diffuse(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    std::size_t nCoeffs = fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> tmp{outarray.size()};
    for (std::size_t i = 0; i < outarray.size(); ++i)
    {
        tmp[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};
    }

    DiffusionLDGET::DiffuseCoeffs(nConvectiveFields, fields, inarray, tmp, pFwd,
                                  pBwd);

    for (std::size_t i = 0; i < outarray.size(); ++i)
    {
        fields[i]->BwdTrans(tmp[i], outarray[i]);
    }
}

void DiffusionLDGET::v_DiffuseCoeffs(
    const std::size_t nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    if (fields[0]->GetGraph()->GetMovement()->GetMoveFlag()) // i.e. if
                                                             // m_ALESolver
    {
        fields[0]->GetTrace()->GetNormals(m_traceNormals);
    }

    std::size_t nDim      = fields[0]->GetCoordim(0);
    std::size_t nPts      = fields[0]->GetTotPoints();
    std::size_t nCoeffs   = fields[0]->GetNcoeffs();
    std::size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nOutput   = outarray.size();

    TensorOfArray3D<NekDouble> qfield{nDim};
    for (std::size_t j = 0; j < nDim; ++j)
    {
        qfield[j] = Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        for (std::size_t i = 0; i < nConvectiveFields; ++i)
        {
            qfield[j][i] = Array<OneD, NekDouble>{nPts, 0.0};
        }
    }

    DiffuseCalcDerivative(fields, inarray, qfield, pFwd, pBwd);

    // Initialize viscous tensor
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> viscTensor{nDim};
    for (std::size_t j = 0; j < nDim; ++j)
    {
        viscTensor[j] = Array<OneD, Array<OneD, NekDouble>>{nOutput};
        for (std::size_t i = 0; i < nOutput; ++i)
        {
            viscTensor[j][i] = Array<OneD, NekDouble>{nPts, 0.0};
        }
    }

    DiffuseVolumeFlux(fields, inarray, qfield, viscTensor);
    Array<OneD, Array<OneD, NekDouble>> traceflux{nOutput};
    for (std::size_t i = 0; i < nOutput; ++i)
    {
        traceflux[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }
    DiffuseTraceFlux(fields, inarray, qfield, viscTensor, traceflux, pFwd,
                     pBwd);

    Array<OneD, Array<OneD, NekDouble>> qdbase{nDim};
    for (std::size_t i = 0; i < nOutput; ++i)
    {

        for (std::size_t j = 0; j < nDim; ++j)
        {
            qdbase[j] = viscTensor[j][i];
        }

        fields[i]->IProductWRTDerivBase(qdbase, outarray[i]);

        Vmath::Neg(nCoeffs, outarray[i], 1);
        fields[i]->AddTraceIntegral(traceflux[i], outarray[i]);

        fields[i]->SetPhysState(false);
    }

    // If mesh is not distorted
    if (!fields[0]->GetGraph()->GetMovement()->GetMeshDistortedFlag())
    {
        for (std::size_t i = 0; i < outarray.size(); ++i)
        {
            fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
        }
    }
}

void DiffusionLDGET::v_DiffuseCalcDerivative(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    std::size_t nConvectiveFields = fields.size();
    std::size_t nDim              = fields[0]->GetCoordim(0);
    std::size_t nCoeffs           = fields[0]->GetNcoeffs();
    std::size_t nTracePts         = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> tmp{nCoeffs};
    TensorOfArray3D<NekDouble> flux{nDim};
    for (std::size_t j = 0; j < nDim; ++j)
    {
        flux[j] = Array<OneD, Array<OneD, NekDouble>>{nConvectiveFields};
        for (std::size_t i = 0; i < nConvectiveFields; ++i)
        {
            flux[j][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        }
    }

    NumFluxforScalar(fields, inarray, flux, pFwd, pBwd);

    for (std::size_t j = 0; j < nDim; ++j)
    {
        for (std::size_t i = 0; i < nConvectiveFields; ++i)
        {
            fields[i]->IProductWRTDerivBase(j, inarray[i], tmp);
            Vmath::Neg(nCoeffs, tmp, 1);
            fields[i]->AddTraceIntegral(flux[j][i], tmp);
            fields[i]->SetPhysState(false);
            fields[i]->MultiplyByElmtInvMass(tmp, tmp);
            fields[i]->BwdTrans(tmp, qfield[j][i]);
        }
    }
}

void DiffusionLDGET::v_DiffuseVolumeFlux(
    [[maybe_unused]] const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    TensorOfArray3D<NekDouble> &qfield, TensorOfArray3D<NekDouble> &viscTensor,
    [[maybe_unused]] Array<OneD, int> &nonZeroIndex)
{
    m_fluxVector(inarray, qfield, viscTensor);
}

void DiffusionLDGET::v_DiffuseTraceFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] TensorOfArray3D<NekDouble> &qfield,
    TensorOfArray3D<NekDouble> &viscTensor,
    Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd,
    [[maybe_unused]] Array<OneD, int> &nonZeroIndex)
{
    NumFluxforVector(fields, inarray, viscTensor, TraceFlux, pFwd, pBwd);
}

void DiffusionLDGET::NumFluxforScalar(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &ufield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &uflux,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    std::size_t nTracePts  = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nvariables = fields.size();
    std::size_t nDim       = fields[0]->GetCoordim(0);

    Array<OneD, NekDouble> fluxtemp{nTracePts, 0.0};
    Array<OneD, Array<OneD, NekDouble>> numflux{nvariables};
    // Get the sign of (v \cdot n), v = an arbitrary vector
    // Evaluate upwind flux:
    // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
    for (std::size_t i = 0; i < nvariables; ++i)
    {
        // Upwind
        numflux[i] = {pFwd[i]};

        // Imposing weak boundary condition with flux
        if (fields[0]->GetBndCondExpansions().size())
        {
            ApplyScalarBCs(fields, i, ufield[i], pFwd[i], pBwd[i], numflux[i]);
        }

        for (std::size_t j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nTracePts, m_traceNormals[j], 1, numflux[i], 1,
                        uflux[j][i], 1);
        }
    }
}

void DiffusionLDGET::ApplyScalarBCs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const std::size_t var,
    [[maybe_unused]] const Array<OneD, const NekDouble> &ufield,
    const Array<OneD, const NekDouble> &Fwd,
    [[maybe_unused]] const Array<OneD, const NekDouble> &Bwd,
    Array<OneD, NekDouble> &penaltyflux)
{
    // Number of boundary regions
    std::size_t nBndRegions = fields[var]->GetBndCondExpansions().size();
    std::size_t cnt         = 0;
    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Number of boundary expansion related to that region
        std::size_t nBndEdges =
            fields[var]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts = fields[var]
                                          ->GetBndCondExpansions()[i]
                                          ->GetExp(e)
                                          ->GetTotPoints();

            std::size_t id1 =
                fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);

            std::size_t id2 = fields[0]->GetTrace()->GetPhys_Offset(
                fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Wall") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Symmetry") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallViscous") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallAdiabatic") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallRotational"))
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }
            // For Dirichlet boundary condition: uflux = g_D
            else if (fields[var]
                         ->GetBndConditions()[i]
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Vmath::Vcopy(
                    nBndEdgePts,
                    &(fields[var]->GetBndCondExpansions()[i]->GetPhys())[id1],
                    1, &penaltyflux[id2], 1);
            }
            // For Neumann boundary condition: uflux = u+
            else if ((fields[var]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }
        }
    }
}

/**
 * @brief Build the numerical flux for the 2nd order derivatives
 * todo: add variable coeff and h dependence to penalty term
 */
void DiffusionLDGET::NumFluxforVector(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &ufield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, NekDouble>> &qflux,
    const Array<OneD, Array<OneD, NekDouble>> &uFwd,
    const Array<OneD, Array<OneD, NekDouble>> &uBwd)
{
    std::size_t nTracePts  = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nvariables = qflux.size();
    std::size_t nDim       = qfield.size();

    Array<OneD, NekDouble> qFwd{nTracePts};
    Array<OneD, NekDouble> qBwd{nTracePts};
    Array<OneD, NekDouble> qfluxtemp{nTracePts, 0.0};
    Array<OneD, NekDouble> uterm{nTracePts};

    Array<OneD, Array<OneD, NekDouble>> fluxPen{nvariables};
    for (std::size_t i = 0; i < nvariables; ++i)
    {
        fluxPen[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    // Evaulate upwind flux:
    // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
    for (std::size_t i = 0; i < nvariables; ++i)
    {
        Vmath::Vsub(nTracePts, uFwd[i], 1, uBwd[i], 1, uterm, 1);
        Vmath::Smul(nTracePts, -m_C11, uterm, 1, uterm, 1);
        qflux[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        for (std::size_t j = 0; j < nDim; ++j)
        {
            //  Compute Fwd and Bwd value of ufield of jth direction
            fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);

            // Downwind
            Vmath::Vcopy(nTracePts, qBwd, 1, qfluxtemp, 1);

            Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                        qfluxtemp, 1);

            // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
            Vmath::Vadd(nTracePts, fluxPen[i], 1, qfluxtemp, 1, qfluxtemp, 1);

            // Imposing weak boundary condition with flux
            if (fields[0]->GetBndCondExpansions().size())
            {
                ApplyVectorBCs(fields, i, j, qfield[j][i], qFwd, qBwd,
                               qfluxtemp);
            }

            // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
            // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
            // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
            Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1, qflux[i], 1);
        }
    }
}

/**
 * Diffusion: Imposing weak boundary condition for q with flux
 *  uflux = g_D  on Dirichlet boundary condition
 *  uflux = u_Fwd  on Neumann boundary condition
 */
void DiffusionLDGET::ApplyVectorBCs(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const std::size_t var, const std::size_t dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &qfield,
    const Array<OneD, const NekDouble> &qFwd,
    [[maybe_unused]] const Array<OneD, const NekDouble> &qBwd,
    Array<OneD, NekDouble> &penaltyflux)
{
    std::size_t nBndRegions = fields[var]->GetBndCondExpansions().size();
    std::size_t cnt         = 0;

    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }
        std::size_t nBndEdges =
            fields[var]->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts = fields[var]
                                          ->GetBndCondExpansions()[i]
                                          ->GetExp(e)
                                          ->GetTotPoints();

            std::size_t id1 =
                fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);

            std::size_t id2 = fields[0]->GetTrace()->GetPhys_Offset(
                fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Wall") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "Symmetry") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallViscous") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallAdiabatic") ||
                boost::iequals(
                    fields[var]->GetBndConditions()[i]->GetUserDefined(),
                    "WallRotational"))
            {
                Vmath::Zero(nBndEdgePts, &penaltyflux[id2], 1);
            }
            // For Dirichlet boundary condition:
            // qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
            else if (fields[var]
                         ->GetBndConditions()[i]
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Vmath::Vmul(nBndEdgePts, &m_traceNormals[dir][id2], 1,
                            &qFwd[id2], 1, &penaltyflux[id2], 1);
            }
            // For Neumann boundary condition: qflux = g_N
            else if ((fields[var]->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vmul(
                    nBndEdgePts, &m_traceNormals[dir][id2], 1,
                    &(fields[var]->GetBndCondExpansions()[i]->GetPhys())[id1],
                    1, &penaltyflux[id2], 1);
            }
        }
    }
}

} // namespace PENKNIFE
