#pragma once

#include "QuasineutralityKernels.hpp"
#include "QuasineutralityOp.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename TData>
class QuasineutralityOpImpl : public QuasineutralityOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;

public:
    QuasineutralityOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                          const std::vector<std::string> &components)
        : QuasineutralityOp<TData>(expansionList, components)
    {
    }

    static std::string className;

    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<QuasineutralityOpImpl<ExecSpace, TData>>(
            expansionList, components);
    }

protected:
    void v_SetSpecies(const std::vector<unsigned int> &densityIds,
                      const std::vector<unsigned int> &momentumIds,
                      const std::vector<TData> &charges,
                      const std::vector<TData> &masses,
                      TData electronMass) override
    {
        m_nions = densityIds.size();
        // Precompute the ion momentum weights once, outside Apply.
        m_electronMass = electronMass;
        std::vector<TData> momentumWeights(m_nions);
        for (unsigned int s = 0; s < m_nions; ++s)
        {
            momentumWeights[s] = (electronMass / masses[s]) * charges[s];
            ASSERTL0(std::isfinite(momentumWeights[s]),
                     "Quasineutrality momentum weights must be finite.");
        }
        if (m_nions != 0)
        {
            m_densityIds =
                LibUtilities::MemoryRegion<unsigned int>::template FromVector<
                    MemSpace, unsigned int>(densityIds);
            m_momentumIds =
                LibUtilities::MemoryRegion<unsigned int>::template FromVector<
                    MemSpace, unsigned int>(momentumIds);
            m_momentumWeights = LibUtilities::MemoryRegion<
                TData>::template FromVector<MemSpace, TData>(momentumWeights);
            m_charges = LibUtilities::MemoryRegion<TData>::template FromVector<
                MemSpace, TData>(charges);
        }
    }

    void v_Apply(
        LibUtilities::Field<TData, FieldState::Phys> &in,
        LibUtilities::Field<TData, FieldState::Phys> &parallelCurrent,
        LibUtilities::Field<TData, FieldState::Phys> &electronDensity,
        LibUtilities::Field<TData, FieldState::Phys> &electronMomentum) override
    {
        // This is pointwise: carry the input's interleave layout to the output.
        // Validate before changing the output's layout or obtaining pointers.
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const auto &input    = in.GetBlocks()[blk];
            const auto &output   = electronDensity.GetBlocks()[blk];
            const auto &momentum = electronMomentum.GetBlocks()[blk];
            const auto &current  = parallelCurrent.GetBlocks()[blk];
            ASSERTL0(input.CompSize() == output.CompSize() &&
                         input.GetNumData() == output.GetNumData() &&
                         input.GetNumElements() == output.GetNumElements() &&
                         input.CompSize() == momentum.CompSize() &&
                         input.GetNumData() == momentum.GetNumData() &&
                         input.GetNumElements() == momentum.GetNumElements() &&
                         input.CompSize() == current.CompSize() &&
                         input.GetNumData() == current.GetNumData() &&
                         input.GetNumElements() == current.GetNumElements() &&
                         input.GetInterleaveWidth() ==
                             current.GetInterleaveWidth(),
                     "Quasineutrality input and output blocks differ.");
        }
        electronDensity.SetInterleaveWidth(in);
        electronMomentum.SetInterleaveWidth(in);

        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            auto &input             = in.GetBlocks()[blk];
            auto &output            = electronDensity.GetBlocks()[blk];
            const auto *densityIds =
                m_nions == 0 ? nullptr
                             : m_densityIds.template GetPtr<MemSpace, ReadOnly>(
                                   streamID);
            auto &momentum = electronMomentum.GetBlocks()[blk];
            const auto *momentumIds =
                m_nions == 0
                    ? nullptr
                    : m_momentumIds.template GetPtr<MemSpace, ReadOnly>(
                          streamID);
            const auto *momentumWeights =
                m_nions == 0
                    ? nullptr
                    : m_momentumWeights.template GetPtr<MemSpace, ReadOnly>(
                          streamID);
            const auto *charges =
                m_nions == 0
                    ? nullptr
                    : m_charges.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto &current       = parallelCurrent.GetBlocks()[blk];
            const size_t stride = input.CompSize() * in.GetNumHomoModes();
            QuasineutralKernel<ExecSpace>(
                stride, m_nions, densityIds, momentumIds, charges,
                momentumWeights, QuasineutralityOp<TData>::InvalidId, stride,
                stride, m_electronMass,
                input.template GetPtr<MemSpace, ReadOnly>(streamID),
                current.template GetPtr<MemSpace, ReadOnly>(streamID),
                output.template GetPtr<MemSpace, WriteOnly>(streamID),
                momentum.template GetPtr<MemSpace, WriteOnly>(streamID),
                streamID);
        }
    }

private:
    unsigned int m_nions = 0;
    TData m_electronMass = 0;
    LibUtilities::MemoryRegion<unsigned int> m_densityIds;
    LibUtilities::MemoryRegion<TData> m_charges;
    LibUtilities::MemoryRegion<unsigned int> m_momentumIds;
    LibUtilities::MemoryRegion<TData> m_momentumWeights;
};
} // namespace PENKNIFE::detail
