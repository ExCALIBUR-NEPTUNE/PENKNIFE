#include "DoubleDiffusiveField.hpp"
namespace PENKNIFE
{
std::string DoubleDiffusiveField::class_name =
    SC::GetEquationSystemFactory().RegisterCreatorFunction(
        "DoubleDiffusiveField", DoubleDiffusiveField::create,
        "Anisotropic density and heat diffusion");
DoubleDiffusiveField::DoubleDiffusiveField(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : SingleDiffusiveField(session, graph, false)
{
}
} // namespace PENKNIFE
