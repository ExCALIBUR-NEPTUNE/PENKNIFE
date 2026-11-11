#pragma once
#include "PlasmaSystem.hpp"
namespace PENKNIFE
{
class NullSystem : public PlasmaSystem
{
public:
    friend class MemoryManager<NullSystem>;
    static std::string class_name;
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        auto p = MemoryManager<NullSystem>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }
    ~NullSystem() override = default;

protected:
    NullSystem(const LU::SessionReaderSharedPtr &session,
               const SD::MeshGraphSharedPtr &graph)
        : PlasmaSystem(session, graph)
    {
    }
    void v_InitObject(bool declareFields = true) override;
    void DoOdeRhs(PhysField &in, PhysField &out, const double &time,
                  const double &factor);
};
} // namespace PENKNIFE
