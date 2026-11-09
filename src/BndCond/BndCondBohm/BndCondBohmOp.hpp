#pragma once

#include <SolverCore/BndCond/BndCondUpdateOp.hpp>

namespace PENKNIFE
{
using namespace Nektar;

template <typename T> struct BohmSpecies
{
    unsigned count, electronEnergy;
    const unsigned *density, *momentum, *energy;
    const T *mass, *charge;
    T scale;
};

// Apply transforms the seeded interior state into the prescribed physical
// normal flux. These Bohm Neumann entries hold fluxes, not derivatives.
// The plasma trace operator consumes them directly.
template <typename T>
class BndCondBohmOp : public SolverCore::BndCondUpdateOp<T>
{
public:
    static inline const std::string name = "BndCondBohm";
    static std::shared_ptr<BndCondBohmOp> Create(
        const MultiRegions::ExpListSharedPtr &expansion,
        const std::vector<std::string> &components,
        const std::string &exec = "")
    {
        return MultiRegions::Operator<T>::template Create<BndCondBohmOp>(
            expansion, components, exec);
    }
    void SetSpecies(BohmSpecies<T> species)
    {
        m_species = species;
    }
    // Three B components followed by parallel current, each at the boundary
    // storage's component stride. Gathered before the standard Apply call.
    LibUtilities::MemoryRegion<T> &MagneticCurrent(unsigned blk)
    {
        return m_magneticCurrent.at(blk);
    }

protected:
    BohmSpecies<T> m_species{};
    std::vector<LibUtilities::MemoryRegion<T>> m_magneticCurrent;
    BndCondBohmOp(const MultiRegions::ExpListSharedPtr &expansion,
                  const std::vector<std::string> &components)
        : SolverCore::BndCondUpdateOp<T>(expansion, components)
    {
    }
};
} // namespace PENKNIFE
