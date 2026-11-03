#pragma once

#include <MultiRegions/DisContField.h>

#include <LibUtilities/BasicUtils/Field/MemoryRegion.hpp>

#include "BndCondBohmKernels.hpp"
#include "BndCondBohmOp.hpp"

using namespace Nektar;

namespace Nektar::detail
{

/// Tag this operator claims.
inline bool IsPressureOutflowTag(const std::string &tag)
{
    return tag == "PressureOutflow";
}

template <typename ExecSpace, typename EoSParamType, typename TData>
class BndCondBohmOpImpl : public BndCondBohmOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;

public:
    BndCondBohmOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                      const std::vector<std::string> &components)
        : BndCondBohmOp<TData>(expansionList, components)
    {
        m_numComp  = components.size();
        m_coordDim = MultiRegions::GetCollection(expansionList, 0)
                         .GetExpVector()[0]
                         ->GetCoordim();

        ASSERTL0(m_numComp == m_coordDim + 2,
                 "A compressible pressure outflow needs density, one momentum "
                 "component per direction and energy.");

        SetUpEquationOfState(expansionList->GetSession(), m_EoS);
    }

    // className - for OperatorFactory
    static std::string className;

    // instantiation function for CreatorFunction in OperatorFactory
    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<
            BndCondBohmOpImpl<ExecSpace, EoSParamType, TData>>(
            expansionList, components);
    }

protected:
    unsigned int m_numComp  = 0;
    unsigned int m_coordDim = 0;
    EoSParamType m_EoS;

    /// Prescribed static pressure per owned block, taken from the session
    /// before the storage is overwritten. Empty for unowned blocks. Held in
    /// the execution space's memory: the boundary storage lives there, so
    /// both the copy at setup and the reads in the apply kernel stay
    /// on-space.
    std::vector<LibUtilities::MemoryRegion<TData>> m_pressure;

    std::vector<bool> v_GetOwnedBlocks() const override
    {
        ASSERTL1(this->m_bndCondPhysOp,
                 "No boundary storage attached: call SetBndCondPhysOp().");

        const auto numBlocks = this->m_bndCondPhysOp->GetNumBlocks();
        std::vector<bool> owned(numBlocks, false);

        for (unsigned blk = 0; blk < numBlocks; ++blk)
        {
            owned[blk] = IsPressureOutflowTag(
                this->m_bndCondPhysOp->GetBndUserDefined(blk, 0));
        }

        return owned;
    }

    /// Copy the prescribed pressure out of the session-evaluated storage,
    /// while it is still there.
    void v_SetBndCondPhysOp() override
    {
        const auto owned = v_GetOwnedBlocks();
        m_pressure.clear();
        m_pressure.resize(owned.size());

        for (unsigned blk = 0; blk < owned.size(); ++blk)
        {
            if (!owned[blk])
            {
                continue;
            }

            const auto stride  = this->m_bndCondPhysOp->GetBlockCompStride(blk);
            const auto *bndPtr = this->m_bndCondPhysOp->GetBlockPtr(blk);

            // The boundary storage lives in the execution space's memory, so
            // the copy runs there as a kernel - a host-side assign from this
            // pointer is a host dereference of device memory on a real
            // device.
            m_pressure[blk] = LibUtilities::MemoryRegion<TData>(stride);
            m_pressure[blk].template Initialize<MemSpace>(0);
            TData *pOut =
                m_pressure[blk].template GetPtr<MemSpace, ReadWrite>();

            CopyPrescribedPressureBlock<ExecSpace>(stride, m_numComp - 1,
                                                   bndPtr, pOut);
        }
    }

    void v_Apply() override
    {
        ASSERTL1(this->m_bndCondPhysOp,
                 "No boundary storage attached: call SetBndCondPhysOp().");

        const auto owned = v_GetOwnedBlocks();

        for (unsigned blk = 0; blk < owned.size(); ++blk)
        {
            if (!owned[blk])
            {
                continue;
            }

            // Stream 0, as the trace flux work around this runs there.
            const unsigned int streamID = 0;

            const auto stride = this->m_bndCondPhysOp->GetBlockCompStride(blk);
            auto bndPtr = this->m_bndCondPhysOp->UpdateBlockPtr(blk, streamID);
            const TData *pOut =
                m_pressure[blk].template GetPtr<MemSpace, ReadOnly>(streamID);
            const TData *norms = this->m_bndCondPhysOp->GetBlockNormals(blk);

            ASSERTL1(m_pressure[blk].size() == stride,
                     "Prescribed outflow pressure was taken for a different "
                     "block size than the storage now reports.");

            ApplyPressureOutflowBlock<ExecSpace>(stride, m_coordDim,
                                                 m_numComp - 1, m_EoS, pOut,
                                                 norms, bndPtr, streamID);
        }
    }
};

} // namespace Nektar::detail
