#include "NullSystem.hpp"
namespace PENKNIFE
{
std::string NullSystem::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "NullSystem", NullSystem::create,
        "Particles in a prescribed plasma background");
void NullSystem::v_InitObject(bool declareFields)
{
    PlasmaSystem::v_InitObject(declareFields);
    m_timeOp->DefineExplicitRhs(&NullSystem::DoOdeRhs, this);
}
void NullSystem::DoOdeRhs(PhysField &in, PhysField &out, const double &time,
                          const double &factor)
{
    m_math.zero(out);
    // The fluid background is prescribed; reaction feedback is diagnostic only.
    for (auto &forcing : m_forcing)
    {
        forcing->SetAppend(true);
        forcing->Apply(in, out, time, factor);
    }
}
} // namespace PENKNIFE
