#pragma once

#include "FluxVolAdvKernels.hpp"
#include "FluxVolAdvOp.hpp"

namespace PENKNIFE::detail
{

template <typename ExecSpace, typename TData>
class FluxVolAdvOpImpl : public FluxVolAdvOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;

public:
    FluxAdvVolOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                     const std::vector<std::string> &components)
        : FluxVolAdvOp<TData>(std::move(expansionList), components)
    {
        m_nDim = expansionList->GetCoordim(0);
    }

    // className - for OperatorFactory
    static std::string className;

    // instantiation function for CreatorFunction in OperatorFactory
    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<FluxVolAdvOpImpl<ExecSpace, TData>>(
            expansionList, components);
    }

protected:
    unsigned int m_nDim;

    void v_Apply(LibUtilities::Field<TData, FieldState::Phys> &in,
                 LibUtilities::Field<TData, FieldState::Phys> &out) override
    {
        for (unsigned int blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const unsigned int streamID = blk + 1;

            auto &inblock  = in.GetBlocks()[blk];
            auto &outblock = out.GetBlocks()[blk];
            auto &velblock = this->m_advVel->GetBlocks()[blk];
            auto &Bblock   = this->m_b_unit.GetBlocks()[blk];

            auto inbase = inblock.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto outbase =
                outblock.template GetPtr<MemSpace, WriteOnly>(streamID);
            auto velbase =
                velblock.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto Bbase = Bblock.template GetPtr<MemSpace, ReadOnly>(streamID);

            auto n_ids = m_n_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto v_ids = m_v_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto e_ids = m_e_ids.template GetPtr<MemSpace, ReadOnly>(streamID);
            auto masses =
                m_masses.template GetPtr<MemSpace, ReadOnly>(streamID);

            const auto npts      = outblock.CompSize();
            const auto inStride  = inblock.CompSize();
            const auto outStride = outblock.CompSize();
            const auto velStride = velblock.CompSize();
            const size_t BStride = Bblock.CompSize();

            FluxVolAdvKernel<ExecSpace>(npts, m_n_species, n_ids, v_ids, e_ids,
                                        masses, m_nDim, velStride, inStride,
                                        outStride, BStride, velbase, inbase,
                                        outbase, Bbase, streamID);
        }
    }
};

} // namespace PENKNIFE::detail