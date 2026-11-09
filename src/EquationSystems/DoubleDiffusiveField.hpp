#pragma once
#include "SingleDiffusiveField.hpp"
namespace PENKNIFE
{
class DoubleDiffusiveField : public SingleDiffusiveField
{
public:
    friend class MemoryManager<DoubleDiffusiveField>;
    static std::string class_name;
    static SC::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        auto p = MemoryManager<DoubleDiffusiveField>::AllocateSharedPtr(session,
                                                                        graph);
        p->InitObject();
        return p;
    }
    ~DoubleDiffusiveField() override = default;

protected:
    DoubleDiffusiveField(const LU::SessionReaderSharedPtr &session,
                         const SD::MeshGraphSharedPtr &graph);
    std::string DiffusionModel() const override
    {
        return "DoubleDiffusive";
    }
};
} // namespace PENKNIFE
