#include "BndCond/BndCondBohm/BndCondBohmOp.hpp"
#include "FluxTraceAdvDiv/FluxTraceAdvDiffOp.hpp"
#include <MultiRegions/DisContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <iostream>
#include <stdexcept>
using namespace Nektar;
using Field = LibUtilities::Field<double, FieldState::Phys>;
void fill(Field &f, const std::vector<double> &v)
{
    f.Initialize<NektarSpaces::HostSpace>(0);
    for (auto &block : f.GetBlocks())
    {
        auto *p = block.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (unsigned c = 0; c < v.size(); ++c)
            std::fill(p + c * block.CompSize(), p + (c + 1) * block.CompSize(),
                      v[c]);
    }
}
void close(double a, double b)
{
    if (!std::isfinite(a) ||
        std::abs(a - b) > 1e-8 * std::max(1., std::abs(b)))
        throw std::runtime_error("Mismatch " + std::to_string(a) + " vs " +
                                 std::to_string(b));
}
int main(int argc, char **argv)
{
    auto session = LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph = SpatialDomains::MeshGraphIO::Read(session);
    auto dg = MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(
        session, graph, "n", true, true, Collections::eNoCollection);
    dg->SetDataWarehouse();
    auto attrs =
        MultiRegions::GetLocTraceBlockAttributes<double, FieldState::Phys>(dg);
    Field state("state", attrs, 4, 1), gradient("gradient", attrs, 8, 1),
        b("b", attrs, 3, 1), current("j", attrs, 1, 1),
        out("out", attrs, 4, 1), nf("nf", attrs, 4, 1);
    fill(state, {2, 0, 3, 3});
    fill(gradient, {0, 0, 0, 0, 0, 0, 0, 0});
    fill(b, {1, 0, 0});
    fill(current, {0});
    fill(nf, {0, 0, 0, 0});
    auto op = PENKNIFE::FluxTraceAdvDiffOp<double>::Create(
        dg, session->GetVariables(), "Serial");
    op->SetIndices({0}, {1}, {2}, {1}, {1});
    op->SetElectronEnergyIndex(3);
    op->SetMagneticTrace(b);
    op->SetParallelCurrentTrace(current);
    op->AddBndCondUpdateOp(PENKNIFE::BndCondBohmOp<double>::Create(
        dg, session->GetVariables(), "Serial"));
    std::vector<double> baseline;
    for (unsigned form = 0; form < 2; ++form)
    {
        if (form)
        {
            op->SetViscousTraceForm(
                SolverCore::TraceFluxOp<
                    double>::ViscousTraceForm::AverageOfNormalFluxes);
            op->SetBoundaryGradientTrace(&gradient);
        }
        op->UpdateBndCond(state);
        auto &second = form ? nf : gradient;
        op->BeginParallelExchange(state, 0);
        op->BeginParallelExchange(second, 1);
        op->Apply(state, second, out);
        op->EndParallelExchange();
        op->ApplyParallel(state, second, out);
        auto a = out.ToArray();
        size_t stride = a.size() / 4;
        if (!form)
            baseline.assign(a.begin(), a.end());
        else
            for (size_t i = 0; i < a.size(); ++i)
                close(a[i], baseline[i]);
        for (unsigned c : {0u, 2u, 3u})
        {
            double loss = c == 0 ? 4. : c == 2 ? 10. : 18.;
            double maxLoss = 0;
            for (size_t i = 0; i < stride; ++i)
            {
                close(std::isfinite(a[c * stride + i]), 1.);
                double x = a[c * stride + i];
                close(x * (x + loss),
                      0.); // zero interior/tangent or outward Bohm advection
                maxLoss = std::max(maxLoss, -x);
            }
            session->GetComm()->AllReduce(maxLoss, LibUtilities::ReduceMax);
            close(maxLoss, loss);
        }
    }
    std::cout
        << "Bohm operator integration passed: actual tagged boundaries and "
           "both trace forms.\n";
    session->Finalise();
}
