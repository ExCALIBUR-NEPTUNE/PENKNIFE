#pragma once

#include "Misc/Constants.hpp"
#include <LibUtilities/BasicUtils/Field/MemoryRegion.hpp>
#include <cmath>
#include <limits>
#include <vector>

namespace PENKNIFE
{
template <typename TData, typename MemSpace> class SpeciesOperator
{
public:
    static constexpr unsigned InvalidId = std::numeric_limits<unsigned>::max();

    void SetSpecies(const std::vector<unsigned> &density,
                      const std::vector<unsigned> &momentum,
                      const std::vector<unsigned> &energy,
                      const std::vector<TData> &masses,
                      const std::vector<TData> &charges)
    {
        ASSERTL0(density.size() == momentum.size() &&
                     density.size() == energy.size() &&
                     density.size() == masses.size() &&
                     density.size() == charges.size(),
                 "Species metadata vectors must be aligned.");

        n_species = density.size();

        if (n_species)
        {
            n_ids = Nektar::LibUtilities::MemoryRegion<
                unsigned>::template FromVector<MemSpace, unsigned>(density);
            v_ids = Nektar::LibUtilities::MemoryRegion<
                unsigned>::template FromVector<MemSpace, unsigned>(momentum);
            e_ids = Nektar::LibUtilities::MemoryRegion<
                unsigned>::template FromVector<MemSpace, unsigned>(energy);
            m_masses = Nektar::LibUtilities::MemoryRegion<
                TData>::template FromVector<MemSpace, TData>(masses);
            m_charges = Nektar::LibUtilities::MemoryRegion<
                TData>::template FromVector<MemSpace, TData>(charges);
        }
        m_species_configured = true;
    }
    void SetElectronEnergyId(unsigned id)
    {
        ee_id = id;
    }

protected:
    unsigned n_species        = 0;
    bool m_species_configured = false;
    std::vector<bool> m_speciesComponents;
    Nektar::LibUtilities::MemoryRegion<unsigned> n_ids, v_ids, e_ids;
    Nektar::LibUtilities::MemoryRegion<TData> m_masses, m_charges;

    unsigned ee_id        = InvalidId;
    TData m_electron_mass = constants::m_e_m_p;
};
} // namespace PENKNIFE
