#pragma once

#include <MultiRegions/Common/Operator.hpp>
#include <cmath>
#include <limits>

namespace PENKNIFE
{
using namespace Nektar;

// Quasineutral electron quantities derived from the current plasma state.
// Parallel current determines electron momentum alongside quasineutral density.
template <typename TData>
class QuasineutralityOp : public MultiRegions::Operator<TData>
{
public:
    static std::shared_ptr<QuasineutralityOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        return MultiRegions::Operator<TData>::template Create<
            QuasineutralityOp>(expansionList, components, execStr);
    }

    static inline const std::string name = "Quasineutrality";

    static constexpr unsigned int InvalidId =
        std::numeric_limits<unsigned int>::max();

    // Compact ion-only metadata. InvalidId denotes no evolved ion momentum;
    // that species contributes density but zero momentum. Masses share the
    // plasma normalization (proton mass units).
    // The implementation copies the metadata into its execution memory space.
    void SetSpecies(const std::vector<unsigned int> &densityIds,
                    const std::vector<unsigned int> &momentumIds,
                    const std::vector<TData> &charges,
                    const std::vector<TData> &masses, TData electronMass)
    {
        ASSERTL0(densityIds.size() == charges.size() &&
                     densityIds.size() == momentumIds.size() &&
                     densityIds.size() == masses.size(),
                 "Quasineutrality needs aligned ion density, momentum, charge "
                 "and mass metadata.");
        ASSERTL0(std::isfinite(electronMass) && electronMass > 0,
                 "Quasineutrality electron mass must be finite and positive.");
        for (size_t s = 0; s < densityIds.size(); ++s)
        {
            ASSERTL0(densityIds[s] < this->m_components.size(),
                     "Quasineutrality ion density index is outside the input "
                     "state.");
            ASSERTL0(momentumIds[s] == InvalidId ||
                         momentumIds[s] < this->m_components.size(),
                     "Quasineutrality ion momentum index is outside the input "
                     "state.");
            ASSERTL0(std::isfinite(masses[s]) && masses[s] > 0,
                     "Quasineutrality ion masses must be finite and positive.");
            ASSERTL0(std::isfinite(charges[s]),
                     "Quasineutrality ion charges must be finite.");
        }
        v_SetSpecies(densityIds, momentumIds, charges, masses, electronMass);
        m_speciesConfigured = true;
    }

    // in and parallelCurrent are read-only; both scalar outputs are
    // overwritten. M_e = m_e (sum_i Z_i M_i/m_i - j_parallel), with M_i = m_i
    // n_i u_i. Current is normalized by e Nnorm cs: positive current reduces
    // M_e. Its block interleave layout must match the supplied stage field.
    void Apply(LibUtilities::Field<TData, FieldState::Phys> &in,
               LibUtilities::Field<TData, FieldState::Phys> &parallelCurrent,
               LibUtilities::Field<TData, FieldState::Phys> &electronDensity,
               LibUtilities::Field<TData, FieldState::Phys> &electronMomentum)
    {
        ASSERTL0(m_speciesConfigured,
                 "Call QuasineutralityOp::SetSpecies before Apply.");
        ASSERTL0(&in != &electronDensity && &in != &electronMomentum &&
                     &electronDensity != &electronMomentum &&
                     &parallelCurrent != &electronDensity &&
                     &parallelCurrent != &electronMomentum,
                 "Quasineutrality requires separate inputs, density and "
                 "momentum fields.");
        ASSERTL0(
            in.GetNumComponents() == this->m_components.size() &&
                parallelCurrent.GetNumComponents() == 1 &&
                electronDensity.GetNumComponents() == 1 &&
                electronMomentum.GetNumComponents() == 1,
            "Quasineutrality needs the plasma state, scalar current, density "
            "and momentum outputs.");
        ASSERTL0(
            in.GetNumHomoModes() == parallelCurrent.GetNumHomoModes() &&
                in.GetBlocks().size() == parallelCurrent.GetBlocks().size() &&
                in.GetNumHomoModes() == electronDensity.GetNumHomoModes() &&
                in.GetBlocks().size() == electronDensity.GetBlocks().size() &&
                in.GetNumHomoModes() == electronMomentum.GetNumHomoModes() &&
                in.GetBlocks().size() == electronMomentum.GetBlocks().size(),
            "Quasineutrality input and output field layouts differ.");
        v_Apply(in, parallelCurrent, electronDensity, electronMomentum);
    }

    void operator()(
        LibUtilities::Field<TData, FieldState::Phys> &in,
        LibUtilities::Field<TData, FieldState::Phys> &parallelCurrent,
        LibUtilities::Field<TData, FieldState::Phys> &electronDensity,
        LibUtilities::Field<TData, FieldState::Phys> &electronMomentum)
    {
        Apply(in, parallelCurrent, electronDensity, electronMomentum);
    }

protected:
    QuasineutralityOp(const MultiRegions::ExpListSharedPtr &expansionList,
                      const std::vector<std::string> &components)
        : MultiRegions::Operator<TData>(expansionList, components)
    {
    }
    ~QuasineutralityOp() override = default;

    virtual void v_SetSpecies(const std::vector<unsigned int> &densityIds,
                              const std::vector<unsigned int> &momentumIds,
                              const std::vector<TData> &charges,
                              const std::vector<TData> &masses,
                              TData electronMass) = 0;
    virtual void v_Apply(
        LibUtilities::Field<TData, FieldState::Phys> &in,
        LibUtilities::Field<TData, FieldState::Phys> &parallelCurrent,
        LibUtilities::Field<TData, FieldState::Phys> &electronDensity,
        LibUtilities::Field<TData, FieldState::Phys> &electronMomentum) = 0;

private:
    bool m_speciesConfigured = false;
};
} // namespace PENKNIFE
