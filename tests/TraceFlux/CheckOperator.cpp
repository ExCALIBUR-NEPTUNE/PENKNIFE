#include "FluxTraceAdvDiv/FluxTraceAdvDiffOp.hpp"
#include <MultiRegions/DisContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <iostream>
#include <stdexcept>
using namespace Nektar;
using Field = LibUtilities::Field<double, FieldState::Phys>;
void close(double a, double b)
{
    if (!std::isfinite(a) || std::abs(a - b) > 1e-9 * std::max(1., std::abs(b)))
        throw std::runtime_error("Mismatch " + std::to_string(a) +
                                 " expected " + std::to_string(b));
}
void fill(Field &f, const std::vector<double> &values)
{
    f.Initialize<NektarSpaces::HostSpace>(0);
    for (auto &block : f.GetBlocks())
    {
        double *p = block.GetPtr<NektarSpaces::HostSpace, ReadWrite>();
        for (unsigned c = 0; c < values.size(); ++c)
            std::fill(p + c * block.CompSize(), p + (c + 1) * block.CompSize(),
                      values[c]);
    }
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
        normalFlux("normalFlux", attrs, 4, 1), b("b", attrs, 3, 1),
        current("current", attrs, 1, 1), out("out", attrs, 4, 1);
    fill(state, {2, 4, 3, 3});
    fill(gradient, {0, 0, 0, 0, 0, 0, 0, 0});
    fill(normalFlux, {0, 0, 0, 0});
    fill(b, {1, 0, 0});
    fill(current, {.3});
    auto op = PENKNIFE::FluxTraceAdvDiffOp<double>::Create(
        dg, session->GetVariables(), "Serial");
    op->SetIndices({0}, {1}, {2}, {1}, {1});
    op->SetElectronEnergyIndex(3);
    op->SetMagneticTrace(b);
    op->SetParallelCurrentTrace(current);
    auto apply = [&](Field &second)
    {
        op->BeginParallelExchange(state, 0);
        op->BeginParallelExchange(second, 1);
        op->Apply(state, second, out);
        op->EndParallelExchange();
        op->ApplyParallel(state, second, out);
    };
    apply(gradient);
    auto baseline = out.ToArray();
    for (unsigned c = 0; c < 4; ++c)
    {
        double maxValue = 0, sum = 0;
        const auto stride = baseline.size() / 4;
        for (size_t i = 0; i < stride; ++i)
        {
            close(std::isfinite(baseline[c * stride + i]), 1);
            maxValue = std::max(maxValue, std::abs(baseline[c * stride + i]));
            sum += baseline[c * stride + i];
        }
        // Every nonzero entry is a physical x-face flux; y-face fluxes are
        // zero.
        close(maxValue, std::vector<double>{4, 10, 6, 5.55}[c]);
        close(sum, 0);
    }
    op->SetViscousTraceForm(SolverCore::TraceFluxOp<
                            double>::ViscousTraceForm::AverageOfNormalFluxes);
    op->SetBoundaryGradientTrace(&gradient);
    apply(normalFlux);
    auto averaged = out.ToArray();
    for (size_t i = 0; i < baseline.size(); ++i)
        close(averaged[i], baseline[i]);
    op->SetAppend(true);
    apply(normalFlux);
    auto appended = out.ToArray();
    for (size_t i = 0; i < baseline.size(); ++i)
        close(appended[i], 2 * baseline[i]);
    // Current is sampled again each cycle, even when the field pointer is
    // unchanged.
    op->SetAppend(false);
    fill(current, {2.3});
    apply(normalFlux);
    auto changed = out.ToArray();
    const auto stride = changed.size() / 4;
    for (size_t i = 0; i < stride; ++i)
    {
        for (unsigned c = 0; c < 3; ++c)
            close(changed[c * stride + i], baseline[c * stride + i]);
        close(changed[3 * stride + i],
              baseline[3 * stride + i] * (2.55 / 5.55));
    }
    if (session->GetComm()->GetSize() > 1)
    {
        const double rank = session->GetComm()->GetRank();
        fill(state, {2 + rank, 4 * (1 + rank), 3 * (1 + rank), 3 * (1 + rank)});
        fill(current, {.3 + .4 * rank});
        apply(normalFlux);
        auto discontinuous = out.ToArray();
        const size_t span = discontinuous.size() / 4;
        for (unsigned c = 0; c < 4; ++c)
        {
            double sum = 0, magnitude = 0;
            for (size_t i = 0; i < span; ++i)
            {
                sum += discontinuous[c * span + i];
                magnitude += std::abs(discontinuous[c * span + i]);
            }
            session->GetComm()->AllReduce(sum, LibUtilities::ReduceSum);
            session->GetComm()->AllReduce(magnitude, LibUtilities::ReduceSum);
            close(sum / std::max(1., magnitude), 0);
        }
    }
    if (session->GetComm()->GetRank() == 0)
        std::cout << "Trace operator integration passed: both forms, "
                     "conservation, append and refreshed current.\n";
    session->Finalise();
}
