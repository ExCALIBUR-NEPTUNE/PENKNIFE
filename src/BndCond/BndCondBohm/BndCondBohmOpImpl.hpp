#pragma once

#include "BndCondBohmKernels.hpp"

namespace PENKNIFE::detail
{
template <typename ExecSpace, typename T>
class BndCondBohmOpImpl : public BndCondBohmOp<T>
{
    using MemSpace = typename ExecSpace::memory_space;

public:
    BndCondBohmOpImpl(const MultiRegions::ExpListSharedPtr &expansion,
                      const std::vector<std::string> &components)
        : BndCondBohmOp<T>(expansion, components), m_comps(components.size()),
          m_dim(expansion->GetCoordim(0))
    {
        const auto session = expansion->GetSession();
        session->LoadParameter("BohmAcousticGamma", m_parameters.acousticGamma,
                               3.);
        session->LoadParameter("SheathIonTransmission",
                               m_parameters.ionTransmission, 2.5);
        session->LoadParameter("SheathElectronTransmission",
                               m_parameters.electronTransmission, 4.5);
        for (T value :
             {m_parameters.acousticGamma, m_parameters.ionTransmission,
              m_parameters.electronTransmission})
            ASSERTL0(std::isfinite(value) && value > 0,
                     "Bohm coefficients must be finite and positive.");
    }
    static std::string className;
    static std::unique_ptr<MultiRegions::Operator<T>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansion,
        const std::vector<std::string> &components)
    {
        return std::make_unique<BndCondBohmOpImpl>(expansion, components);
    }

protected:
    unsigned m_comps, m_dim;
    BohmParameters<T> m_parameters;
    std::vector<bool> v_GetOwnedBlocks() const override
    {
        const auto blocks = this->m_bndCondPhysOp->GetNumBlocks();
        std::vector<bool> owned(blocks, false);
        for (unsigned blk = 0; blk < blocks; ++blk)
            owned[blk] =
                this->m_bndCondPhysOp->GetBndUserDefined(blk, 0) == "Bohm";
        return owned;
    }
    void v_SetBndCondPhysOp() override
    {
        const auto owned = v_GetOwnedBlocks();
        this->m_magneticCurrent.resize(owned.size());
        for (unsigned blk = 0; blk < owned.size(); ++blk)
        {
            if (!owned[blk])
                continue;
            for (unsigned c = 0; c < m_comps; ++c)
                ASSERTL0(this->m_bndCondPhysOp->GetBndUserDefined(blk, c) ==
                                 "Bohm" &&
                             this->m_bndCondPhysOp->GetBndCondType(blk, c) ==
                                 SpatialDomains::eNeumann,
                         "Tag every component of a Bohm region as Neumann with "
                         "USERDEFINEDTYPE=Bohm.");
            const auto stride = this->m_bndCondPhysOp->GetBlockCompStride(blk);
            this->m_magneticCurrent[blk] =
                LibUtilities::MemoryRegion<T>(4 * stride);
            this->m_magneticCurrent[blk].template Initialize<MemSpace>(0);
        }
    }
    T v_BndFluxWeight(unsigned, unsigned) const override
    {
        return 0;
    }
    void v_Apply() override
    {
        ASSERTL0(this->m_species.count &&
                     this->m_species.electronEnergy < m_comps,
                 "Set Bohm species metadata before Apply.");
        const auto owned = v_GetOwnedBlocks();
        for (unsigned blk = 0; blk < owned.size(); ++blk)
        {
            if (!owned[blk])
                continue;
            const auto stride = this->m_bndCondPhysOp->GetBlockCompStride(blk);
            ApplyBohmBoundary<ExecSpace>(
                this->m_species, m_parameters, stride, m_dim,
                this->m_bndCondPhysOp->GetBlockNormals(blk),
                this->m_magneticCurrent[blk]
                    .template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                this->m_bndCondPhysOp->UpdateBlockPtr(blk, this->m_opStreamID),
                this->m_opStreamID);
        }
    }
};
} // namespace PENKNIFE::detail
