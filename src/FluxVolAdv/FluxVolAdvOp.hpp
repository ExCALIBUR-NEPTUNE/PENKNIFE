#pragma once

#include "SolverCore/Advection/AdvectionVolumeFluxOp.hpp"

namespace PENKNIFE
{

// Upwind base class
// Defines the apply operator to enforce apply parameter types
template <typename TData>
class FluxVolAdvOp : public SolverCore::AdvectionVolumeFluxOp<TData>
{
public:
    static std::shared_ptr<FluxVolAdvOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        return std::dynamic_pointer_cast<FluxVolAdvOp<TData>>(
            SolverCore::AdvectionVolumeFluxOp<TData>::Create(
                expansionList, components, name, execStr));
    }

    static inline const std::string name = "FluxVolAdv";

    /// The magnetic field is referenced, not adopted, so that it stays usable
    /// by whoever owns it; it must outlive this operator.
    void SetMagneticField(LibUtilities::Field<TData, FieldState::Phys> &b_unit)
    {
        this->m_b_unit = &b_unit;
    }
    SetIndices(std::vector<unsigned int> &n_ids,
               std::vector<unsigned int> &v_ids,
               std::vector<unsigned int> &e_ids, std::vector<TData> &masses,
               std::vector<TData> &charges)
    {
        m_n_species = n_ids.size();
        m_n_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(n_ids);
        m_v_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(v_ids);
        m_e_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(e_ids);
        m_masses = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(masses);
        m_charges = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(charges);
    }

protected:
    unsigned int m_n_species;
    LibUtilities::MemoryRegion<unsigned int> m_n_ids;
    LibUtilities::MemoryRegion<unsigned int> m_v_ids;
    LibUtilities::MemoryRegion<unsigned int> m_e_ids;
    LibUtilities::MemoryRegion<TData> m_masses;
    LibUtilities::MemoryRegion<TData> m_charges;
    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit = nullptr;

    FluxVolAdvOp(const MultiRegions::ExpListSharedPtr &expansionList,
                 const std::vector<std::string> &components)
        : SolverCore::AdvectionVolumeFluxOp<TData>(expansionList, components)
    {
    }

    ~FluxVolAdvOp() override = default;
};

} // namespace PENKNIFE
