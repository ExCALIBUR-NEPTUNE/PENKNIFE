#pragma once

#include "FluxVolAdvDiffKernels.hpp"
#include "FluxVolAdvDiffOp.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename ClosureType, typename TData>
class FluxVolAdvDiffOpImpl : public FluxVolAdvDiffOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;

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
    void v_SetIndices(const std::vector<unsigned> &n_ids,
                      const std::vector<unsigned> &v_ids,
                      const std::vector<unsigned> &e_ids,
                      const std::vector<TData> &masses,
                      const std::vector<TData> &charges) override
    {
        m_nSpecies = n_ids.size();
        if (m_nSpecies == 0)
            return;
        m_nIds = LibUtilities::MemoryRegion<unsigned>::template FromVector<
            MemSpace, unsigned>(n_ids);
        m_vIds = LibUtilities::MemoryRegion<unsigned>::template FromVector<
            MemSpace, unsigned>(v_ids);
        m_eIds = LibUtilities::MemoryRegion<unsigned>::template FromVector<
            MemSpace, unsigned>(e_ids);
        m_masses = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(masses);
        m_charges = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(charges);
    }
    void v_Apply(
        LibUtilities::Field<TData, FieldState::Phys> &in,
        LibUtilities::Field<TData, FieldState::Phys> &deriv,
        LibUtilities::Field<TData, FieldState::Phys> &out,
        LibUtilities::Field<TData, FieldState::Phys> &totalOut) override
    {
        ASSERTL0(this->m_indicesConfigured && this->m_b_unit && this->m_ne &&
                     this->m_Me && this->m_electronEnergyId != this->InvalidId,
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
            if (m_nSpecies != 0)
            {
                nIds   = m_nIds.template GetPtr<MemSpace, ReadOnly>(streamID);
                vIds   = m_vIds.template GetPtr<MemSpace, ReadOnly>(streamID);
                eIds   = m_eIds.template GetPtr<MemSpace, ReadOnly>(streamID);
                masses = m_masses.template GetPtr<MemSpace, ReadOnly>(streamID);
                charges =
                    m_charges.template GetPtr<MemSpace, ReadOnly>(streamID);
            }
            // Treat all homogeneous planes of a component as a single span.
            const size_t stride = ib.CompSize() * in.GetNumHomoModes();
            FluxVolAdvDiffKernel<ExecSpace, ClosureType>(
                m_parameters, stride, m_dimension, in.GetNumComponents(),
                m_nSpecies, this->m_electronEnergyId, nIds, vIds, eIds, masses,
                charges, this->m_electronMass, m_anomD, m_anomNu, m_anomChi,
                stride, stride, stride, stride, stride,
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

private:
    unsigned m_dimension, m_nSpecies = 0;
    typename ClosureType::Parameters m_parameters;
    TData m_anomD, m_anomNu, m_anomChi;
    LibUtilities::MemoryRegion<unsigned> m_nIds, m_vIds, m_eIds;
    LibUtilities::MemoryRegion<TData> m_masses, m_charges;
};
} // namespace PENKNIFE::detail
