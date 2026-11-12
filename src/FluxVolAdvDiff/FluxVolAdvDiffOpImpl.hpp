#pragma once

#include "../SpeciesOperator.hpp"
#include "FluxVolAdvDiffKernels.hpp"
#include "FluxVolAdvDiffOp.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename ClosureType, typename TData>
class FluxVolAdvDiffOpImpl
    : public FluxVolAdvDiffOp<TData>,
      public SpeciesOperator<TData, typename ExecSpace::memory_space>
{
    using MemSpace = typename ExecSpace::memory_space;

    using SpeciesBase = SpeciesOperator<TData, MemSpace>;
    using SpeciesBase::e_ids;
    using SpeciesBase::ee_id;
    using SpeciesBase::InvalidId;
    using SpeciesBase::m_charges;
    using SpeciesBase::m_masses;
    using SpeciesBase::n_ids;
    using SpeciesBase::n_species;
    using SpeciesBase::v_ids;

public:
    FluxVolAdvDiffOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                         const std::vector<std::string> &components)
        : FluxVolAdvDiffOp<TData>(expansionList, components)
    {
        m_dimension        = expansionList->GetCoordim(0);
        const auto session = expansionList->GetSession();
        session->LoadParameter("Nnorm", m_parameters.Nnorm, 1e18);
        session->LoadParameter("Tnorm", m_parameters.Tnorm, 100.);
        session->LoadParameter("mesh_length", m_parameters.mesh_length, 1.);
        TData Bnorm;
        session->LoadParameter("Bnorm", Bnorm, 1.);
        m_parameters.omega_c = TData(constants::qeomp) * Bnorm;
        for (auto value : {m_parameters.Nnorm, m_parameters.Tnorm,
                           m_parameters.mesh_length, m_parameters.omega_c})
            ASSERTL0(std::isfinite(value) && value > 0,
                     "Braginskii normalization parameters must be finite and "
                     "positive.");
        session->LoadParameter("Anomalous_D", m_anomD, 0.);
        session->LoadParameter("Anomalous_nu", m_anomNu, 0.);
        session->LoadParameter("Anomalous_chi", m_anomChi, 0.);
    }
    static std::string className;
    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<FluxVolAdvDiffOpImpl>(expansionList,
                                                      components);
    }

protected:
    void v_Apply(
        LibUtilities::Field<TData, FieldState::Phys> &in,
        LibUtilities::Field<TData, FieldState::Phys> &deriv,
        LibUtilities::Field<TData, FieldState::Phys> &out,
        LibUtilities::Field<TData, FieldState::Phys> &totalOut) override
    {
        ASSERTL0(this->m_species_configured && this->m_b_unit && this->m_ne &&
                     this->m_Me && this->ee_id != this->InvalidId,
                 "Configure combined-flux species, unit B, electron energy, "
                 "density and momentum before Apply.");
        auto &b  = *this->m_b_unit;
        auto &ne = *this->m_ne;
        auto &Me = *this->m_Me;
        ASSERTL0(&out != &totalOut,
                 "Viscous and total outputs must be separate.");
        for (auto *output : {&out, &totalOut})
            for (auto *input : {&in, &deriv, &b, &ne, &Me})
                ASSERTL0(output != input,
                         "Combined flux output must not alias an input.");
        ASSERTL0(in.GetNumComponents() == this->m_components.size() &&
                     deriv.GetNumComponents() ==
                         in.GetNumComponents() * m_dimension &&
                     out.GetNumComponents() == deriv.GetNumComponents() &&
                     b.GetNumComponents() == 3 && ne.GetNumComponents() == 1 &&
                     Me.GetNumComponents() == 1 &&
                     totalOut.GetNumComponents() == out.GetNumComponents(),
                 "Diffusion needs state, component-major gradients/fluxes and "
                 "three-component unit B.");
        for (auto *field : {&deriv, &out, &totalOut, &b, &ne, &Me})
            ASSERTL0(field->GetNumHomoModes() == in.GetNumHomoModes() &&
                         field->GetBlocks().size() == in.GetBlocks().size(),
                     "Diffusion field block layouts differ.");
        // Validate all blocks before touching the output.
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const auto &input = in.GetBlocks()[blk];
            for (auto *field : {&deriv, &out, &totalOut, &b, &ne, &Me})
            {
                const auto &block = field->GetBlocks()[blk];
                ASSERTL0(block.CompSize() == input.CompSize() &&
                             block.GetNumData() == input.GetNumData() &&
                             block.GetNumElements() == input.GetNumElements(),
                         "Diffusion field blocks differ.");
            }
            for (auto *field : {&deriv, &b, &ne, &Me})
                ASSERTL0(field->GetBlocks()[blk].GetInterleaveWidth() ==
                             input.GetInterleaveWidth(),
                         "Combined flux inputs must use the same interleave "
                         "layout.");
        }
        // The combined operator overwrites both outputs, as required by the
        // AdvDiffVolumeFluxOp interface (no append mode).
        out.SetInterleaveWidth(in);
        totalOut.SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            auto &ib                = in.GetBlocks()[blk];
            auto &db                = deriv.GetBlocks()[blk];
            auto &ob                = out.GetBlocks()[blk];
            auto &bb                = b.GetBlocks()[blk];
            auto &tb                = totalOut.GetBlocks()[blk];
            auto &neb               = ne.GetBlocks()[blk];
            auto &Meb               = Me.GetBlocks()[blk];
            const unsigned *nIds = nullptr, *vIds = nullptr, *eIds = nullptr;
            const TData *masses = nullptr, *charges = nullptr;
            if (this->n_species != 0)
            {
                nIds =
                    this->n_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
                vIds =
                    this->v_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
                eIds =
                    this->e_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
                masses = this->m_masses.template GetPtr<MemSpace, ReadOnly>(
                    streamID);
                charges = this->m_charges.template GetPtr<MemSpace, ReadOnly>(
                    streamID);
            }
            // Treat all homogeneous planes of a component as a single span.
            const size_t stride = ib.CompSize() * in.GetNumHomoModes();
            FluxVolAdvDiffKernel<ExecSpace, ClosureType>(
                m_parameters, stride, m_dimension, in.GetNumComponents(),
                this->n_species, this->ee_id, nIds, vIds, eIds, masses, charges,
                this->m_electron_mass, m_anomD, m_anomNu, m_anomChi, stride,
                stride, stride, stride, stride,
                ib.template GetPtr<MemSpace, ReadOnly>(streamID),
                db.template GetPtr<MemSpace, ReadOnly>(streamID),
                ob.template GetPtr<MemSpace, WriteOnly>(streamID),
                tb.template GetPtr<MemSpace, WriteOnly>(streamID),
                bb.template GetPtr<MemSpace, ReadOnly>(streamID),
                neb.template GetPtr<MemSpace, ReadOnly>(streamID),
                Meb.template GetPtr<MemSpace, ReadOnly>(streamID),
                this->m_scale, streamID);
        }
    }

    void v_SetSpecies(const std::vector<unsigned> &density,
                      const std::vector<unsigned> &momentum,
                      const std::vector<unsigned> &energy,
                      const std::vector<TData> &masses,
                      const std::vector<TData> &charges) override
    {
        SpeciesOperator<TData, MemSpace>::SetSpecies(density, momentum, energy,
                                                     masses, charges);
    }

    void v_SetElectronEnergyId(unsigned id) override
    {
        SpeciesOperator<TData, MemSpace>::SetElectronEnergyId(id);
    }

private:
    unsigned m_dimension;
    typename ClosureType::Parameters m_parameters;
    TData m_anomD, m_anomNu, m_anomChi;
};
} // namespace PENKNIFE::detail
