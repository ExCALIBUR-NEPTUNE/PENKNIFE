#pragma once
#include "PlasmaPointwiseOp.hpp"
#include <LibUtilities/BasicUtils/Utils/UtilsKernels.hpp>
#include <LibUtilities/LoopExecution/LoopExecution.hpp>
#include <cmath>
namespace PENKNIFE::detail
{
using namespace Nektar;
namespace LU = Nektar::LibUtilities;
// Kernels index scalar DoFs directly for every interleave width, including
// AVX storage. No temporary flattening or vector-width conversion is needed.
template <typename ExecSpace, typename TData>
class PlasmaPointwiseOpImpl : public PlasmaPointwiseOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;
    using Field    = typename PlasmaPointwiseOp<TData>::Field;
    static constexpr unsigned InvalidId = PlasmaPointwiseOp<TData>::InvalidId;

public:
    PlasmaPointwiseOpImpl(const MultiRegions::ExpListSharedPtr &expansion,
                          const std::vector<std::string> &components)
        : PlasmaPointwiseOp<TData>(expansion, components),
          m_coordim(expansion->GetCoordim(0))
    {
    }
    static std::string className;
    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansion,
        const std::vector<std::string> &components)
    {
        return std::make_unique<PlasmaPointwiseOpImpl>(expansion, components);
    }
    void SetSpecies(const std::vector<unsigned> &density,
                    const std::vector<unsigned> &momentum,
                    const std::vector<unsigned> &energy,
                    const std::vector<TData> &masses,
                    const std::vector<TData> &charges,
                    const std::vector<unsigned> &ionIds,
                    unsigned electronEnergy, TData electronMass,
                    TData backgroundTemperature) override
    {
        ASSERTL0(density.size() == momentum.size() &&
                     density.size() == energy.size() &&
                     density.size() == masses.size() &&
                     density.size() == charges.size(),
                 "Pointwise species metadata must be aligned.");
        ASSERTL0(std::isfinite(electronMass) && electronMass > 0,
                 "Pointwise electron mass must be finite and positive.");
        for (const auto ion : ionIds)
        {
            ASSERTL0(ion < density.size() && density[ion] != InvalidId &&
                         density[ion] < this->m_components.size(),
                     "Quasineutrality ion density index is outside the input "
                     "state.");
            ASSERTL0(momentum[ion] == InvalidId ||
                         momentum[ion] < this->m_components.size(),
                     "Quasineutrality ion momentum index is outside the input "
                     "state.");
            ASSERTL0(std::isfinite(masses[ion]) && masses[ion] > 0,
                     "Quasineutrality ion masses must be finite and positive.");
            ASSERTL0(std::isfinite(charges[ion]),
                     "Quasineutrality ion charges must be finite.");
        }
        n_species               = density.size();
        n_ions                  = ionIds.size();
        ee_idx                  = electronEnergy;
        me                      = electronMass;
        m_backgroundTemperature = backgroundTemperature;
        if (n_species)
        {
            n_ids = LU::MemoryRegion<unsigned>::template FromVector<MemSpace,
                                                                    unsigned>(
                density);
            v_ids = LU::MemoryRegion<unsigned>::template FromVector<MemSpace,
                                                                    unsigned>(
                momentum);
            e_ids = LU::MemoryRegion<unsigned>::template FromVector<MemSpace,
                                                                    unsigned>(
                energy);
            m_masses =
                LU::MemoryRegion<TData>::template FromVector<MemSpace, TData>(
                    masses);
            m_charges =
                LU::MemoryRegion<TData>::template FromVector<MemSpace, TData>(
                    charges);
        }
        if (n_ions)
            ions = LU::MemoryRegion<unsigned>::template FromVector<MemSpace,
                                                                   unsigned>(
                ionIds);
    }
    void SetSources(const std::vector<unsigned> &density,
                    const std::vector<unsigned> &momentum,
                    const std::vector<unsigned> &energy,
                    unsigned electronEnergy) override
    {
        m_sourceElectron = electronEnergy;
        if (n_species)
        {
            m_sourceDensity = LU::MemoryRegion<unsigned>::template FromVector<
                MemSpace, unsigned>(density);
            m_sourceMomentum = LU::MemoryRegion<unsigned>::template FromVector<
                MemSpace, unsigned>(momentum);
            m_sourceEnergy = LU::MemoryRegion<unsigned>::template FromVector<
                MemSpace, unsigned>(energy);
        }
    }
    void Quasineutrality(Field &in, Field &parallelCurrent, Field &ne,
                         Field &Me) override
    {
        parallelCurrent.SetInterleaveWidth(in);
        ne.SetInterleaveWidth(in);
        Me.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *ions =
                (this->n_ions == 0
                     ? nullptr
                     : this->ions.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *charges =
                (this->n_species == 0
                     ? nullptr
                     : this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_ions = this->n_ions;
            const auto me     = this->me;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto current = parallelCurrent.GetBlocks()[blk]
                               .template GetPtr<MemSpace, ReadOnly>(streamID);
            auto density =
                ne.GetBlocks()[blk].template GetPtr<MemSpace, WriteOnly>(
                    streamID);
            auto momentum =
                Me.GetBlocks()[blk].template GetPtr<MemSpace, WriteOnly>(
                    streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    TData electronDensity  = 0;
                    TData electronMomentum = -me * current[p];
                    for (unsigned ion = 0; ion < n_ions; ++ion)
                    {
                        const unsigned s = ions[ion];
                        electronDensity +=
                            charges[s] * u[n_ids[s] * stride + p];
                        if (v_ids[s] != InvalidId)
                            electronMomentum += (me / masses[s]) * charges[s] *
                                                u[v_ids[s] * stride + p];
                    }
                    density[p]  = electronDensity;
                    momentum[p] = electronMomentum;
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void ParticleFields(Field &in, Field &ne, Field &Me, Field &b_unit,
                        Field &B, Field &E, Field &m_particleFields) override
    {
        ne.SetInterleaveWidth(in);
        Me.SetInterleaveWidth(in);
        b_unit.SetInterleaveWidth(in);
        B.SetInterleaveWidth(in);
        E.SetInterleaveWidth(in);
        m_particleFields.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *e_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->e_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_species               = this->n_species;
            const auto ee_idx                  = this->ee_idx;
            const auto me                      = this->me;
            const auto m_backgroundTemperature = this->m_backgroundTemperature;

            auto &block         = in.GetBlocks()[blk];
            const size_t stride = block.CompSize() * in.GetNumHomoModes();
            auto u = block.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto n = ne.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto momentum =
                Me.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto b =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto magnetic =
                B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto electric =
                E.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto values = m_particleFields.GetBlocks()[blk]
                              .template GetPtr<MemSpace, WriteOnly>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    values[p] = n[p];
                    values[stride + p] =
                        n[p] > 0 && ee_idx != InvalidId
                            ? (2. / 3.) * u[ee_idx * stride + p] / n[p]
                            : m_backgroundTemperature;
                    for (unsigned d = 0; d < 3; ++d)
                    {
                        values[(2 + d) * stride + p] =
                            n[p] > 0
                                ? momentum[p] * b[d * stride + p] / (me * n[p])
                                : 0;
                        values[(5 + d) * stride + p] = electric[d * stride + p];
                        values[(8 + d) * stride + p] = magnetic[d * stride + p];
                    }
                    for (unsigned s = 0; s < n_species; ++s)
                    {
                        const size_t c         = 11 + 5 * s;
                        const TData ns         = (n_ids[s] != InvalidId)
                                                     ? u[n_ids[s] * stride + p]
                                                     : 0;
                        values[c * stride + p] = ns;
                        values[(c + 1) * stride + p] =
                            ns > 0 && (e_ids[s] != InvalidId)
                                ? (2. / 3.) * u[e_ids[s] * stride + p] / ns
                                : m_backgroundTemperature;
                        for (unsigned d = 0; d < 3; ++d)
                            values[(c + 2 + d) * stride + p] =
                                ns > 0 && (v_ids[s] != InvalidId)
                                    ? u[v_ids[s] * stride + p] *
                                          b[d * stride + p] / (m_masses[s] * ns)
                                    : 0;
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void AddSources(Field &out, Field &m_particleSources, Field &b_unit,
                    TData factor, TData timestep) override
    {
        m_particleSources.SetInterleaveWidth(out);
        b_unit.SetInterleaveWidth(out);
        for (unsigned blk = 0; blk < out.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *e_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->e_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_sourceDensity =
                (this->n_species == 0
                     ? nullptr
                     : this->m_sourceDensity
                           .template GetPtr<MemSpace, ReadOnly>(streamID));
            const auto *m_sourceMomentum =
                (this->n_species == 0
                     ? nullptr
                     : this->m_sourceMomentum
                           .template GetPtr<MemSpace, ReadOnly>(streamID));
            const auto *m_sourceEnergy =
                (this->n_species == 0
                     ? nullptr
                     : this->m_sourceEnergy.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_species        = this->n_species;
            const auto ee_idx           = this->ee_idx;
            const auto m_sourceElectron = this->m_sourceElectron;

            const auto stride =
                out.GetBlocks()[blk].CompSize() * out.GetNumHomoModes();
            auto source = m_particleSources.GetBlocks()[blk]
                              .template GetPtr<MemSpace, ReadOnly>(streamID);
            auto b =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto rhs =
                out.GetBlocks()[blk].template GetPtr<MemSpace, ReadWrite>(
                    streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned s = 0; s < n_species; ++s)
                    {
                        if ((n_ids[s] != InvalidId))
                            rhs[n_ids[s] * stride + p] +=
                                factor *
                                source[m_sourceDensity[s] * stride + p] /
                                timestep;
                        if ((e_ids[s] != InvalidId))
                            rhs[e_ids[s] * stride + p] +=
                                factor *
                                source[m_sourceEnergy[s] * stride + p] /
                                timestep;
                        if ((v_ids[s] != InvalidId))
                            for (unsigned d = 0; d < 3; ++d)
                                rhs[v_ids[s] * stride + p] +=
                                    factor * b[d * stride + p] *
                                    source[(m_sourceMomentum[s] + d) * stride +
                                           p] /
                                    timestep;
                    }
                    if (m_sourceElectron != InvalidId)
                        rhs[ee_idx * stride + p] +=
                            factor * source[m_sourceElectron * stride + p] /
                            timestep;
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void ParallelVelocities(Field &in, Field &b_unit, Field &ne, Field &Me,
                            Field &m_velocity) override
    {
        b_unit.SetInterleaveWidth(in);
        ne.SetInterleaveWidth(in);
        Me.SetInterleaveWidth(in);
        m_velocity.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_species = this->n_species;
            const auto me        = this->me;
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto b =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto n = ne.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto M = Me.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto v = m_velocity.GetBlocks()[blk]
                         .template GetPtr<MemSpace, WriteOnly>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned s = 0; s <= n_species; ++s)
                    {
                        TData velocity = 0;
                        if (s == n_species)
                            velocity = n[p] > 0 ? M[p] / (me * n[p]) : 0;
                        else if ((n_ids[s] != InvalidId) &&
                                 (v_ids[s] != InvalidId) &&
                                 u[n_ids[s] * stride + p] > 0)
                            velocity = u[v_ids[s] * stride + p] /
                                       (m_masses[s] * u[n_ids[s] * stride + p]);
                        for (unsigned d = 0; d < m_coordim; ++d)
                            v[(s * m_coordim + d) * stride + p] =
                                velocity * b[d * stride + p];
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void PressureWork(Field &in, Field &out, Field &m_velocityGradient,
                      Field &E, Field &b_unit, TData factor) override
    {
        out.SetInterleaveWidth(in);
        m_velocityGradient.SetInterleaveWidth(in);
        E.SetInterleaveWidth(in);
        b_unit.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *e_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->e_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_charges =
                (this->n_species == 0
                     ? nullptr
                     : this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_species = this->n_species;
            const auto ee_idx    = this->ee_idx;
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto gradient = m_velocityGradient.GetBlocks()[blk]
                                .template GetPtr<MemSpace, ReadOnly>(streamID);
            auto rhs =
                out.GetBlocks()[blk].template GetPtr<MemSpace, ReadWrite>(
                    streamID);
            auto electric =
                E.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto unit =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            // The pressure gradient is in the conservative momentum flux.
            // The remaining parallel Lorentz force is +Z n E_parallel.
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    TData parallelElectric = 0;
                    for (unsigned d = 0; d < 3; ++d)
                        parallelElectric +=
                            unit[d * stride + p] * electric[d * stride + p];
                    for (unsigned s = 0; s < n_species; ++s)
                        if ((n_ids[s] != InvalidId) && (v_ids[s] != InvalidId))
                            rhs[v_ids[s] * stride + p] +=
                                factor * m_charges[s] *
                                u[n_ids[s] * stride + p] * parallelElectric;
                    for (unsigned s = 0; s <= n_species; ++s)
                    {
                        if (s < n_species && (e_ids[s] == InvalidId))
                            continue;
                        const unsigned energy =
                            s == n_species ? ee_idx : e_ids[s];
                        TData divergence = 0;
                        for (unsigned d = 0; d < m_coordim; ++d)
                            divergence +=
                                gradient[((s * m_coordim + d) * m_coordim + d) *
                                             stride +
                                         p];
                        rhs[energy * stride + p] -= factor * (2. / 3.) *
                                                    u[energy * stride + p] *
                                                    divergence;
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void PotentialForcing(Field &in, Field &b_unit, Field &mag_B,
                          Field &m_phiRhs, Field &m_perpCoeff,
                          unsigned omega_idx) override
    {
        b_unit.SetInterleaveWidth(in);
        mag_B.SetInterleaveWidth(in);
        m_phiRhs.SetInterleaveWidth(in);
        m_perpCoeff.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *ions =
                (this->n_ions == 0
                     ? nullptr
                     : this->ions.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_ions    = this->n_ions;
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto b =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto B2 =
                mag_B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto rhs =
                m_phiRhs.GetBlocks()[blk].template GetPtr<MemSpace, WriteOnly>(
                    streamID);
            auto D = m_perpCoeff.GetBlocks()[blk]
                         .template GetPtr<MemSpace, WriteOnly>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    TData rho = 0;
                    for (unsigned ion = 0; ion < n_ions; ++ion)
                    {
                        const unsigned s = ions[ion];
                        if (n_ids[s] != InvalidId)
                            rho += m_masses[s] * u[n_ids[s] * stride + p];
                    }
                    rhs[p] =
                        rho > 0 ? u[omega_idx * stride + p] * B2[p] / rho : 0;
                    for (unsigned d = 0; d < m_coordim; ++d)
                        for (unsigned j = 0; j <= d; ++j)
                            D[LU::GetDiffCoeffMap(m_coordim,
                                                  d * m_coordim + j) *
                                  stride +
                              p] = (d == j ? 1. : 0.) -
                                   b[d * stride + p] * b[j * stride + p];
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void CorrectPotential(Field &in, Field &phi) override
    {
        phi.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *e_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->e_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_charges =
                (this->n_species == 0
                     ? nullptr
                     : this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *ions =
                (this->n_ions == 0
                     ? nullptr
                     : this->ions.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_ions = this->n_ions;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto potential =
                phi.GetBlocks()[blk].template GetPtr<MemSpace, ReadWrite>(
                    streamID);
            // Retain the existing ion-temperature potential correction.
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    if (n_ions == 0)
                        return;
                    const unsigned s = ions[n_ions - 1];
                    if (n_ids[s] != InvalidId && e_ids[s] != InvalidId &&
                        u[n_ids[s] * stride + p] > 0 && m_charges[s] != 0)
                        potential[p] -=
                            (2. / 3.) * u[e_ids[s] * stride + p] /
                            (m_charges[s] * u[n_ids[s] * stride + p]);
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void ElectricField(Field &m_phiGradient, Field &E) override
    {
        E.SetInterleaveWidth(m_phiGradient);
        for (unsigned blk = 0; blk < E.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                E.GetBlocks()[blk].CompSize() * E.GetNumHomoModes();
            auto gradient = m_phiGradient.GetBlocks()[blk]
                                .template GetPtr<MemSpace, ReadOnly>(streamID);
            auto electric =
                E.GetBlocks()[blk].template GetPtr<MemSpace, WriteOnly>(
                    streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned d = 0; d < 3; ++d)
                        electric[d * stride + p] =
                            d < m_coordim ? -gradient[d * stride + p] : 0;
                });
        }

        Nektar::LoopExecutionSetStreamID(0);
    }
    void DriftVelocities(Field &in, Field &m_gradient, Field &B, Field &mag_B,
                         Field &E, Field &b_unit, Field &m_extraVelocity,
                         Field &m_velocity, unsigned omega_idx) override
    {
        m_gradient.SetInterleaveWidth(in);
        B.SetInterleaveWidth(in);
        mag_B.SetInterleaveWidth(in);
        E.SetInterleaveWidth(in);
        b_unit.SetInterleaveWidth(in);
        m_extraVelocity.SetInterleaveWidth(in);
        m_velocity.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *v_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->v_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *e_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->e_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_charges =
                (this->n_species == 0
                     ? nullptr
                     : this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *ions =
                (this->n_ions == 0
                     ? nullptr
                     : this->ions.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const unsigned n_fields = in.GetNumComponents();
            const auto n_species    = this->n_species;
            const auto n_ions       = this->n_ions;
            const auto ee_idx       = this->ee_idx;
            const auto m_coordim    = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto gradient =
                m_gradient.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto magnetic =
                B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto B2 =
                mag_B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto electric =
                E.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto unit =
                b_unit.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto extra = m_extraVelocity.GetBlocks()[blk]
                             .template GetPtr<MemSpace, WriteOnly>(streamID);
            auto velocity = m_velocity.GetBlocks()[blk]
                                .template GetPtr<MemSpace, WriteOnly>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned f = 0; f < n_fields * m_coordim; ++f)
                        extra[f * stride + p] = 0;
                    for (unsigned f = 0; f < (n_species + 1) * m_coordim; ++f)
                        velocity[f * stride + p] = 0;
                    if (B2[p] <= 0)
                        return;
                    TData exb[3] = {0, 0, 0};
                    for (unsigned d = 0; d < 3; ++d)
                    {
                        const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                        exb[d] = (electric[j * stride + p] *
                                      magnetic[k * stride + p] -
                                  electric[k * stride + p] *
                                      magnetic[j * stride + p]) /
                                 B2[p];
                    }
                    for (unsigned d = 0; d < m_coordim; ++d)
                        extra[(omega_idx * m_coordim + d) * stride + p] =
                            exb[d];
                    for (unsigned s = 0; s <= n_species; ++s)
                    {
                        const bool electron = s == n_species;
                        if (!electron &&
                            ((n_ids[s] == InvalidId) || m_charges[s] == 0))
                            continue;
                        const unsigned energy = electron ? ee_idx : e_ids[s];
                        TData density = 0, parallel = 0;
                        if (electron)
                        {
                            for (unsigned ion = 0; ion < n_ions; ++ion)
                            {
                                const unsigned t = ions[ion];
                                if (n_ids[t] != InvalidId)
                                {
                                    density +=
                                        m_charges[t] * u[n_ids[t] * stride + p];
                                    if ((v_ids[t] != InvalidId))
                                        parallel += m_charges[t] *
                                                    u[v_ids[t] * stride + p] /
                                                    m_masses[t];
                                }
                            }
                            parallel = density > 0 ? parallel / density : 0;
                        }
                        else
                        {
                            density  = u[n_ids[s] * stride + p];
                            parallel = density > 0 && (v_ids[s] != InvalidId)
                                           ? u[v_ids[s] * stride + p] /
                                                 (m_masses[s] * density)
                                           : 0;
                        }
                        if (density <= 0)
                            continue;
                        TData gp[3] = {0, 0, 0};
                        if (energy != InvalidId)
                            for (unsigned d = 0; d < m_coordim; ++d)
                                gp[d] =
                                    (2. / 3.) *
                                    gradient[(energy * m_coordim + d) * stride +
                                             p];
                        const TData charge = electron ? -1. : m_charges[s];
                        for (unsigned d = 0; d < m_coordim; ++d)
                        {
                            const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                            // v_dia = B cross grad(p) / (q n B^2).
                            const TData drift =
                                (magnetic[j * stride + p] * gp[k] -
                                 magnetic[k * stride + p] * gp[j]) /
                                (charge * B2[p] * density);
                            if (!electron)
                            {
                                extra[(n_ids[s] * m_coordim + d) * stride + p] =
                                    exb[d] + drift;
                                if ((v_ids[s] != InvalidId))
                                    extra[(v_ids[s] * m_coordim + d) * stride +
                                          p] = exb[d] + drift;
                                velocity[(s * m_coordim + d) * stride + p] =
                                    parallel * unit[d * stride + p] + exb[d] +
                                    drift;
                            }
                            if (energy != InvalidId)
                                extra[(energy * m_coordim + d) * stride + p] =
                                    exb[d] + (5. / 3.) * drift;
                        }
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void PolarizationTensor(Field &in, Field &m_velocity, Field &B,
                            Field &mag_B, Field &m_tensor,
                            TData omega_c) override
    {
        m_velocity.SetInterleaveWidth(in);
        B.SetInterleaveWidth(in);
        mag_B.SetInterleaveWidth(in);
        m_tensor.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *n_ids =
                (this->n_species == 0
                     ? nullptr
                     : this->n_ids.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *m_masses =
                (this->n_species == 0
                     ? nullptr
                     : this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto *ions =
                (this->n_ions == 0
                     ? nullptr
                     : this->ions.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_ions    = this->n_ions;
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto u = in.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                streamID);
            auto velocity =
                m_velocity.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto magnetic =
                B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto B2 =
                mag_B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto tensor =
                m_tensor.GetBlocks()[blk].template GetPtr<MemSpace, WriteOnly>(
                    streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned c = 0; c < m_coordim * m_coordim; ++c)
                        tensor[c * stride + p] = 0;
                    if (B2[p] <= 0)
                        return;
                    for (unsigned ion = 0; ion < n_ions; ++ion)
                    {
                        const unsigned s = ions[ion];
                        TData v[3]       = {0, 0, 0};
                        for (unsigned d = 0; d < m_coordim; ++d)
                            v[d] = velocity[(s * m_coordim + d) * stride + p];
                        for (unsigned d = 0; d < m_coordim; ++d)
                        {
                            const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                            const TData w = u[n_ids[s] * stride + p] *
                                            (v[j] * magnetic[k * stride + p] -
                                             v[k] * magnetic[j * stride + p]) /
                                            B2[p];
                            for (unsigned t = 0; t < m_coordim; ++t)
                                tensor[(t * m_coordim + d) * stride + p] +=
                                    m_masses[s] / omega_c * w * v[t];
                        }
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void PolarizationFlux(Field &in, Field &m_tensorGradient,
                          Field &m_omegaFlux) override
    {
        m_tensorGradient.SetInterleaveWidth(in);
        m_omegaFlux.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto derivative =
                m_tensorGradient.GetBlocks()[blk]
                    .template GetPtr<MemSpace, ReadOnly>(streamID);
            auto flux = m_omegaFlux.GetBlocks()[blk]
                            .template GetPtr<MemSpace, WriteOnly>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    for (unsigned d = 0; d < m_coordim; ++d)
                    {
                        TData divergence = 0;
                        for (unsigned j = 0; j < m_coordim; ++j)
                            divergence +=
                                derivative[((d * m_coordim + j) * m_coordim +
                                            j) *
                                               stride +
                                           p];
                        flux[d * stride + p] = divergence;
                    }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }
    void CompleteVelocities(Field &in, Field &E, Field &B, Field &mag_B,
                            Field &m_velocity) override
    {
        E.SetInterleaveWidth(in);
        B.SetInterleaveWidth(in);
        mag_B.SetInterleaveWidth(in);
        m_velocity.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            Nektar::LoopExecutionSetStreamID(streamID);
            const auto *m_charges =
                (this->n_species == 0
                     ? nullptr
                     : this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                           streamID));
            const auto n_species = this->n_species;
            const auto m_coordim = this->m_coordim;

            const size_t stride =
                in.GetBlocks()[blk].CompSize() * in.GetNumHomoModes();
            auto electric =
                E.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto magnetic =
                B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto B2 =
                mag_B.GetBlocks()[blk].template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            auto v = m_velocity.GetBlocks()[blk]
                         .template GetPtr<MemSpace, ReadWrite>(streamID);
            Nektar::parallel_for<ExecSpace>(
                0, stride, NEKTAR_LAMBDA(const size_t p) {
                    if (B2[p] > 0)
                        for (unsigned d = 0; d < m_coordim; ++d)
                        {
                            const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                            const TData exb = (electric[j * stride + p] *
                                                   magnetic[k * stride + p] -
                                               electric[k * stride + p] *
                                                   magnetic[j * stride + p]) /
                                              B2[p];
                            for (unsigned s = 0; s <= n_species; ++s)
                                if (s == n_species || m_charges[s] != 0)
                                    v[(s * m_coordim + d) * stride + p] += exb;
                        }
                });
        }
        Nektar::LoopExecutionSetStreamID(0);
    }

private:
    unsigned n_species = 0, n_ions = 0, ee_idx = InvalidId,
             m_sourceElectron = InvalidId;
    unsigned m_coordim;
    TData me = 0, m_backgroundTemperature = 0;
    LU::MemoryRegion<unsigned> n_ids, v_ids, e_ids, ions;
    LU::MemoryRegion<TData> m_masses, m_charges;
    LU::MemoryRegion<unsigned> m_sourceDensity, m_sourceMomentum,
        m_sourceEnergy;
};
} // namespace PENKNIFE::detail
