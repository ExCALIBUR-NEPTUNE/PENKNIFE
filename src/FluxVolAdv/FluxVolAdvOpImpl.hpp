#pragma once

#include "FluxVolAdvKernels.hpp"
#include "FluxVolAdvOp.hpp"
#include "SpeciesOperator.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename TData>
class FluxVolAdvOpImpl
    : public FluxVolAdvOp<TData>,
      public SpeciesOperator<TData, typename ExecSpace::memory_space>
{
    using MemSpace = typename ExecSpace::memory_space;

public:
    FluxVolAdvOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                     const std::vector<std::string> &components)
        : FluxVolAdvOp<TData>(expansionList, components)
    {
        m_dimension = expansionList->GetCoordim(0);
    }
    static std::string className;
    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<FluxVolAdvOpImpl>(expansionList, components);
    }

protected:
    void v_Apply(LibUtilities::Field<TData, FieldState::Phys> &in,
                 LibUtilities::Field<TData, FieldState::Phys> &out) override
    {
        ASSERTL0(this->m_species_configured && this->m_b_unit && this->m_ne &&
                     this->m_Me && this->ee_id != this->InvalidId,
                 "Configure combined-flux species, unit B, electron energy, "
                 "density and momentum before Apply.");
        auto &b  = *this->m_b_unit;
        auto &ne = *this->m_ne;
        auto &Me = *this->m_Me;

        // Validate all blocks before touching the output.
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const auto &input = in.GetBlocks()[blk];
            for (auto *field : {&out, &b, &ne, &Me})
            {
                const auto &block = field->GetBlocks()[blk];
                ASSERTL0(block.CompSize() == input.CompSize() &&
                             block.GetNumData() == input.GetNumData() &&
                             block.GetNumElements() == input.GetNumElements(),
                         "Diffusion field blocks differ.");
            }
            for (auto *field : {&b, &ne, &Me})
                ASSERTL0(field->GetBlocks()[blk].GetInterleaveWidth() ==
                             input.GetInterleaveWidth(),
                         "Combined flux inputs must use the same interleave "
                         "layout.");
        }
        // The combined operator overwrites both outputs, as required by the
        // AdvDiffVolumeFluxOp interface (no append mode).
        out.SetInterleaveWidth(in);

        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            auto &ib                = in.GetBlocks()[blk];
            auto &ob                = out.GetBlocks()[blk];
            auto &bb                = b.GetBlocks()[blk];
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
            FluxVolAdvKernel<ExecSpace>(
                stride, m_dimension, in.GetNumComponents(), this->n_species,
                this->ee_id, nIds, vIds, eIds, masses, charges,
                this->m_electron_mass, stride, stride, stride, stride,
                ib.template GetPtr<MemSpace, ReadOnly>(streamID),
                ob.template GetPtr<MemSpace, WriteOnly>(streamID),
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
};
} // namespace PENKNIFE::detail
