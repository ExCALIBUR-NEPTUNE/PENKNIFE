#include "Closure.hpp"
#include "../EquationSystems/PlasmaSystem.hpp"

namespace PENKNIFE
{
ClosureFactory &GetClosureFactory()
{
    static ClosureFactory instance;
    return instance;
}

Closure::Closure(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim)
    : m_system(pSystem), m_spacedim(spaceDim),
      field_to_index(pSystem.lock()->field_to_index),
      b_unit(pSystem.lock()->b_unit), mag_B(pSystem.lock()->mag_B)
{
    this->n_pts   = pSystem.lock()->n_pts;
    this->Nnorm   = pSystem.lock()->Nnorm;
    this->Tnorm   = pSystem.lock()->Tnorm;
    this->omega_c = pSystem.lock()->omega_c;
}
} // namespace PENKNIFE