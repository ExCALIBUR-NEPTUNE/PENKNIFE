#include "ElectrostaticTurbulence.hpp"
namespace PENKNIFE
{
std::string ElectrostaticTurbulence::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "ElectrostaticTurbulence", ElectrostaticTurbulence::create,
        "Electrostatic turbulence with anisotropic diffusion");
ElectrostaticTurbulence::ElectrostaticTurbulence(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : ReducedBraginskii(session, graph)
{
    omega_idx = m_variables.size();
    m_variables.push_back("w");
    m_fieldTemplates.push_back("w");
    m_nVariables = m_variables.size();
}
void ElectrostaticTurbulence::LoadParams()
{
    PlasmaSystem::LoadParams();
    std::string boussinesq;
    m_session->LoadSolverInfo("Boussinesq Approximation", boussinesq, "On");
    ASSERTL0(boussinesq == "On",
             "Only Boussinesq Approximation=On is currently implemented.");
}
void ElectrostaticTurbulence::v_InitialiseFields()
{
    ReducedBraginskii::v_InitialiseFields();
    m_phiExpansion = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, "phi", false, true);
    m_phiExpansion->SetDataWarehouse();
    auto volume =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_expansionLists[0]);
    auto trace = MR::GetLocTraceBlockAttributes<double, FieldState::Phys>(
        m_expansionLists[0]);
    auto phiAttr =
        MR::GetBlockAttributes<double, FieldState::Phys>(m_phiExpansion);
    phi       = PhysField("phi", phiAttr, 1, 1);
    phi_coeff = CoeffField(
        "phi coefficients",
        MR::GetBlockAttributes<double, FieldState::Coeff>(m_phiExpansion), 1,
        1);
    m_phiOutputRhs = CoeffField(
        "Potential output RHS",
        MR::GetBlockAttributes<double, FieldState::Coeff>(m_phiExpansion), 1,
        1);
    m_phiRhs      = PhysField("Potential forcing", phiAttr, 1, 1);
    m_phiGradient = PhysField("Potential gradient", phiAttr, m_coordim, 1);
    m_perpCoeff   = PhysField("Perpendicular tensor", phiAttr,
                              m_coordim * (m_coordim + 1) / 2, 1);
    m_gradient =
        PhysField("State gradient", volume, m_variables.size() * m_coordim, 1);
    m_extraVelocity      = PhysField("Per-field drift velocities", volume,
                                     m_variables.size() * m_coordim, 1);
    m_extraVelocityTrace = PhysField("Drift velocity traces", trace,
                                     m_variables.size() * m_coordim, 1);
    m_omegaFlux          = PhysField("Vorticity flux", volume, m_coordim, 1);
    m_omegaFluxTrace = PhysField("Vorticity flux traces", trace, m_coordim, 1);
    m_tensor         = PhysField("Polarization velocity tensor", volume,
                                 m_coordim * m_coordim, 1);
    m_tensorGradient = PhysField("Polarization tensor gradient", volume,
                                 m_coordim * m_coordim * m_coordim, 1);
    m_math.zero(phi);
    m_math.zero(phi_coeff);
}
void ElectrostaticTurbulence::v_InitialiseOperators()
{
    ReducedBraginskii::v_InitialiseOperators();
    m_phiSolve  = SC::PoissonSolveOp<double>::Create(m_phiExpansion, {"phi"});
    auto solver = SC::LinearSolverOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiSolve->SetLinearSolver(solver);
    m_phiSolve->SetVarDiffCoeff(m_perpCoeff);
    m_phiSolve->SetPrecon(
        SC::PreconOp<double>::Create(m_phiExpansion, {"phi"}, "Diagonal"));
    m_phiBackward = MR::BwdTransOp<double>::Create(m_phiExpansion, {"phi"});
    m_derivative =
        MR::PhysDerivOp<double>::Create(m_expansionLists[0], m_fieldTemplates);
    m_phiDerivative = MR::PhysDerivOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiOutputInnerProduct =
        MR::IProductWRTBaseOp<double>::Create(m_phiExpansion, {"phi"});
    m_phiOutputInverseMass =
        MR::MultiplyByElmtInvMassOp<double>::Create(m_phiExpansion, {"phi"});
    m_tensorDerivative = MR::PhysDerivOp<double>::Create(
        m_expansionLists[0], m_tensor.GetComponentNames());
    m_advectionFlux->SetExtraFlux(m_extraVelocity, m_omegaFlux, omega_idx);
    m_traceFlux->SetExtraAdvectionTrace(m_extraVelocityTrace, m_omegaFluxTrace,
                                        omega_idx);
}
void ElectrostaticTurbulence::UpdateElectricField(PhysField &in, double time)
{
    for (auto *f : {&b_unit, &mag_B, &m_phiRhs, &m_perpCoeff, &phi})
        f->SetInterleaveWidth(in);
    m_math.zero(m_perpCoeff);
    m_math.zero(m_phiRhs);
    for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
    {
        const size_t stride = in.GetBlocks()[blk].CompSize();
        auto u =
            in.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto b =
            b_unit.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto B2 =
            mag_B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto rhs = m_phiRhs.GetBlocks()[blk]
                       .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        auto D = m_perpCoeff.GetBlocks()[blk]
                     .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (size_t p = 0; p < stride; ++p)
        {
            double rho = 0;
            for (unsigned s : GetIons())
                if (HasDensity(s))
                    rho += m_masses[s] * u[n_ids[s] * stride + p];
            if (rho > 0)
                rhs[p] = u[omega_idx * stride + p] * B2[p] / rho;
            for (unsigned d = 0; d < m_coordim; ++d)
                for (unsigned j = 0; j <= d; ++j)
                    D[LU::GetDiffCoeffMap(m_coordim, d * m_coordim + j) *
                          stride +
                      p] = (d == j ? 1. : 0.) -
                           b[d * stride + p] * b[j * stride + p];
        }
    }
    m_phiSolve->UpdateBndCoeffs(time);
    m_phiSolve->UpdatePrecon();
    m_math.zero(phi_coeff);
    // PoissonSolve uses the physical sign div(D grad(phi)) = rhs.
    m_phiSolve->Apply(m_phiRhs, phi_coeff);
    m_phiBackward->Apply(phi_coeff, phi);
    phi.SetInterleaveWidth(in);
    for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
    {
        const size_t stride = in.GetBlocks()[blk].CompSize();
        auto u =
            in.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto potential =
            phi.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        // Retain the existing ion-temperature potential correction.
        if (!GetIons().empty())
        {
            const unsigned s = GetIons().back();
            if (HasDensity(s) && HasEnergy(s))
                for (size_t p = 0; p < stride; ++p)
                    if (u[n_ids[s] * stride + p] > 0 && m_charges[s] != 0)
                        potential[p] -=
                            (2. / 3.) * u[e_ids[s] * stride + p] /
                            (m_charges[s] * u[n_ids[s] * stride + p]);
        }
    }
    m_phiDerivative->Apply(phi, m_phiGradient);
    E.SetInterleaveWidth(m_phiGradient);
    m_math.zero(E);
    for (unsigned blk = 0; blk < E.GetBlocks().size(); ++blk)
    {
        const size_t stride = E.GetBlocks()[blk].CompSize();
        auto gradient       = m_phiGradient.GetBlocks()[blk]
                            .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto electric =
            E.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (unsigned d = 0; d < m_coordim; ++d)
            for (size_t p = 0; p < stride; ++p)
                electric[d * stride + p] = -gradient[d * stride + p];
    }
}
void ElectrostaticTurbulence::PrepareAdvection(PhysField &in, double time)
{
    ReducedBraginskii::PrepareAdvection(in, time);
    UpdateElectricField(in, time);
    m_derivative->Apply(in, m_gradient);
    for (auto *f : {&B, &mag_B, &E, &b_unit, &m_extraVelocity, &m_velocity})
        f->SetInterleaveWidth(in);
    m_gradient.SetInterleaveWidth(in);
    m_math.zero(m_extraVelocity);
    m_math.zero(m_velocity);
    for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
    {
        const size_t stride = in.GetBlocks()[blk].CompSize();
        auto u =
            in.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto gradient = m_gradient.GetBlocks()[blk]
                            .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto magnetic =
            B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto B2 =
            mag_B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto electric =
            E.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto unit =
            b_unit.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto extra = m_extraVelocity.GetBlocks()[blk]
                         .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        auto velocity = m_velocity.GetBlocks()[blk]
                            .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (size_t p = 0; p < stride; ++p)
        {
            if (B2[p] <= 0)
                continue;
            double exb[3] = {0, 0, 0};
            for (unsigned d = 0; d < 3; ++d)
            {
                const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                exb[d] = (electric[j * stride + p] * magnetic[k * stride + p] -
                          electric[k * stride + p] * magnetic[j * stride + p]) /
                         B2[p];
            }
            for (unsigned d = 0; d < m_coordim; ++d)
                extra[(omega_idx * m_coordim + d) * stride + p] = exb[d];
            for (unsigned s = 0; s <= n_species; ++s)
            {
                const bool electron = s == n_species;
                if (!electron && (!HasDensity(s) || m_charges[s] == 0))
                    continue;
                const unsigned energy = electron ? ee_idx : e_ids[s];
                double density = 0, parallel = 0;
                if (electron)
                {
                    for (unsigned t : GetIons())
                        if (HasDensity(t))
                        {
                            density += m_charges[t] * u[n_ids[t] * stride + p];
                            if (HasMomentum(t))
                                parallel += m_charges[t] *
                                            u[v_ids[t] * stride + p] /
                                            m_masses[t];
                        }
                    parallel = density > 0 ? parallel / density : 0;
                }
                else
                {
                    density = u[n_ids[s] * stride + p];
                    parallel =
                        density > 0 && HasMomentum(s)
                            ? u[v_ids[s] * stride + p] / (m_masses[s] * density)
                            : 0;
                }
                if (density <= 0)
                    continue;
                double gp[3] = {0, 0, 0};
                if (energy != InvalidId)
                    for (unsigned d = 0; d < m_coordim; ++d)
                        gp[d] = (2. / 3.) *
                                gradient[(energy * m_coordim + d) * stride + p];
                const double charge = electron ? -1. : m_charges[s];
                for (unsigned d = 0; d < m_coordim; ++d)
                {
                    const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                    // v_dia = B cross grad(p) / (q n B^2).
                    const double drift = (magnetic[j * stride + p] * gp[k] -
                                          magnetic[k * stride + p] * gp[j]) /
                                         (charge * B2[p] * density);
                    if (!electron)
                    {
                        extra[(n_ids[s] * m_coordim + d) * stride + p] =
                            exb[d] + drift;
                        if (HasMomentum(s))
                            extra[(v_ids[s] * m_coordim + d) * stride + p] =
                                exb[d] + drift;
                        velocity[(s * m_coordim + d) * stride + p] =
                            parallel * unit[d * stride + p] + exb[d] + drift;
                    }
                    if (energy != InvalidId)
                        extra[(energy * m_coordim + d) * stride + p] =
                            exb[d] + (5. / 3.) * drift;
                }
            }
        }
    }
    // omega flux = sum_s m_s/omega_c div[(n_s v_s x B/B^2) tensor v_s].
    m_omegaFlux.SetInterleaveWidth(in);
    m_math.zero(m_omegaFlux);
    for (unsigned s : GetIons())
    {
        m_tensor.SetInterleaveWidth(in);
        m_math.zero(m_tensor);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const size_t stride = in.GetBlocks()[blk].CompSize();
            auto u =
                in.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto velocity = m_velocity.GetBlocks()[blk]
                                .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto magnetic =
                B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto B2 = mag_B.GetBlocks()[blk]
                          .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto tensor = m_tensor.GetBlocks()[blk]
                              .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
            for (size_t p = 0; p < stride; ++p)
            {
                double v[3] = {0, 0, 0};
                for (unsigned d = 0; d < m_coordim; ++d)
                    v[d] = velocity[(s * m_coordim + d) * stride + p];
                if (B2[p] <= 0)
                    continue;
                for (unsigned d = 0; d < m_coordim; ++d)
                {
                    const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                    const double w = u[n_ids[s] * stride + p] *
                                     (v[j] * magnetic[k * stride + p] -
                                      v[k] * magnetic[j * stride + p]) /
                                     B2[p];
                    for (unsigned t = 0; t < m_coordim; ++t)
                        tensor[(t * m_coordim + d) * stride + p] = w * v[t];
                }
            }
        }
        m_tensorDerivative->Apply(m_tensor, m_tensorGradient);
        m_tensorGradient.SetInterleaveWidth(m_omegaFlux);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const size_t stride = in.GetBlocks()[blk].CompSize();
            auto derivative     = m_tensorGradient.GetBlocks()[blk]
                                  .GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto flux = m_omegaFlux.GetBlocks()[blk]
                            .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
            for (size_t p = 0; p < stride; ++p)
                for (unsigned d = 0; d < m_coordim; ++d)
                    for (unsigned j = 0; j < m_coordim; ++j)
                        flux[d * stride + p] +=
                            m_masses[s] / omega_c *
                            derivative[((d * m_coordim + j) * m_coordim + j) *
                                           stride +
                                       p];
        }
    }
    m_traceExtract->Apply(m_extraVelocity, m_extraVelocityTrace);
    m_traceExtract->Apply(m_omegaFlux, m_omegaFluxTrace);
}
void ElectrostaticTurbulence::CompleteVelocities(PhysField &in)
{
    // Pressure work uses ExB plus parallel flow, before diamagnetic heat flux.
    E.SetInterleaveWidth(m_velocity);
    B.SetInterleaveWidth(m_velocity);
    mag_B.SetInterleaveWidth(m_velocity);
    for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
    {
        const size_t stride = in.GetBlocks()[blk].CompSize();
        auto electric =
            E.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto magnetic =
            B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto B2 =
            mag_B.GetBlocks()[blk].GetPtr<NektarSpaces::HostSpace, ReadOnly>();
        auto v = m_velocity.GetBlocks()[blk]
                     .GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (size_t p = 0; p < stride; ++p)
            if (B2[p] > 0)
                for (unsigned s = 0; s <= n_species; ++s)
                    if (s == n_species || m_charges[s] != 0)
                        for (unsigned d = 0; d < m_coordim; ++d)
                        {
                            const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                            v[(s * m_coordim + d) * stride + p] +=
                                (electric[j * stride + p] *
                                     magnetic[k * stride + p] -
                                 electric[k * stride + p] *
                                     magnetic[j * stride + p]) /
                                B2[p];
                        }
    }
}
SC::CoefficientFieldViews ElectrostaticTurbulence::v_GetOutputFieldViews()
{
    auto views = PlasmaSystem::v_GetOutputFieldViews();
    m_phiOutputInnerProduct->Apply(phi, m_phiOutputRhs);
    m_phiOutputInverseMass->Apply(m_phiOutputRhs, phi_coeff);
    views.push_back({phi_coeff, m_phiExpansion, m_coordim, {{"phi", 0}}});
    return views;
}
} // namespace PENKNIFE
