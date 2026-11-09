#include "Closures/DiffusiveFields.hpp"
#include "EquationSystems/DoubleDiffusiveField.hpp"
#include "EquationSystems/ElectrostaticTurbulence.hpp"
#include "EquationSystems/NullSystem.hpp"
#include "EquationSystems/ReducedBraginskii.hpp"
#include "EquationSystems/SingleDiffusiveField.hpp"

#include <SpatialDomains/MeshGraphIO.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>

using namespace Nektar;
using namespace PENKNIFE;
void CheckParticleCoupling(const LibUtilities::SessionReaderSharedPtr &,
                           const SpatialDomains::MeshGraphSharedPtr &,
                           const MultiRegions::ExpListSharedPtr &);
void CheckTransportCoefficients()
{
    using Policy = PENKNIFE::closures::DoubleDiffusive<double>;
    Policy::Parameters parameters{};
    parameters.kpar    = 3;
    parameters.kperp   = 5;
    parameters.kipar   = 7;
    parameters.kiperp  = 11;
    parameters.kicross = 13;
    parameters.kepar   = 17;
    parameters.keperp  = 19;
    parameters.kecross = 23;
    parameters.kci     = 2;
    unsigned density[] = {0}, energy[] = {1};
    double mass[] = {4}, charge[] = {2};
    double state[] = {2, 6, 12}, gradient[] = {1, 0, 6, 0, 6, 0},
           magnetic[] = {1, 0, 1}, flux[6] = {};
    Policy::Evaluate(parameters, 0, 2, 1, 2, density, energy, mass, charge, 1,
                     1, 1, 1, state, gradient, magnetic, flux, 1.);
    const double collisionSum = 8 / std::sqrt(2.);
    const double expectedDensity =
        .5 * (6 * std::pow(2., 2.5) / 8 + 40 / std::sqrt(2.));
    const double expectedHeat =
        .5 * (14 * std::pow(2., 2.5) / (8 * collisionSum) +
              22 * collisionSum / std::sqrt(2.));
    if (std::abs(flux[0] - expectedDensity) > 1e-10 ||
        std::abs(flux[2] - expectedHeat) > 1e-10 ||
        std::abs(flux[3] - 13) > 1e-10 || std::abs(flux[4]) > 1e-10)
        throw std::runtime_error(
            "Density or thermal transport coefficients are incorrect.");
}
template <class Model> class CheckSystem : public Model
{
    using Field = LibUtilities::Field<double, FieldState::Phys>;

public:
    CheckSystem(const LibUtilities::SessionReaderSharedPtr &session,
                const SpatialDomains::MeshGraphSharedPtr &graph)
        : Model(session, graph)
    {
    }
    void Check(bool manufactured)
    {
        this->InitObject();
        this->SetInitialConditions();
        auto attrs = MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
            this->m_expansionLists[0]);
        Field rhs("test rhs", attrs, this->m_variables, 1);
        const double time = 0., factor = .37;
        this->DoOdeRhs(this->m_fields, rhs, time, factor);
        rhs.SetInterleaveWidth(this->m_fields);
        Field coords("coordinates", attrs, 3, 1);
        coords.template Initialize<NektarSpaces::HostSpace>(0.);
        auto coordinateOp = MultiRegions::ExpressionOp<double>::Create(
            this->m_expansionLists[0], {"x", "y", "z"});
        std::vector<LibUtilities::EquationSharedPtr> expressions;
        for (const auto &name : {"x", "y", "z"})
            expressions.push_back(std::make_shared<LibUtilities::Equation>(
                this->m_session->GetInterpreter(), name));
        coordinateOp->SetExpressions(expressions);
        coordinateOp->Apply(coords, coords);
        coords.SetInterleaveWidth(rhs);
        double error = 0., reference = 0., maximum = 0.;
        size_t count           = 0;
        const unsigned density = this->GetDensityId(0);
        for (unsigned b = 0; b < rhs.GetBlocks().size(); ++b)
        {
            auto &block         = rhs.GetBlocks()[b];
            const size_t stride = block.CompSize(),
                         width  = block.GetInterleaveWidth(),
                         np     = block.GetNumData();
            auto actual =
                block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto initial =
                this->m_fields.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto xyz =
                coords.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            for (size_t e = 0; e < block.GetNumElements(); ++e)
                for (size_t q = 0; q < np; ++q)
                {
                    const size_t p =
                        (e / width) * np * width + q * width + e % width;
                    const double expected =
                        manufactured ? -factor * .15 * M_PI * M_PI * .1 *
                                           std::sin(M_PI * xyz[p])
                                     : 0.;
                    const double initialExpected =
                        manufactured ? 2. + .1 * std::sin(M_PI * xyz[p]) : 2.;
                    if (std::abs(initial[density * stride + p] -
                                 initialExpected) > 1e-6)
                        throw std::runtime_error(
                            "Species initial conditions differ: " +
                            std::to_string(initial[density * stride + p]) +
                            " expected " + std::to_string(initialExpected));
                    for (unsigned c = 0; c < this->m_variables.size(); ++c)
                    {
                        const double value = actual[c * stride + p];
                        if (!std::isfinite(value))
                            throw std::runtime_error(
                                "Non-finite equation RHS.");
                        if (!manufactured)
                            maximum = std::max(maximum, std::abs(value));
                    }
                    error +=
                        std::pow(actual[density * stride + p] - expected, 2);
                    reference += expected * expected;
                    ++count;
                }
        }
        if (manufactured && std::sqrt(error / reference) > .015)
            throw std::runtime_error(
                "Manufactured anisotropic diffusion RHS is inaccurate.");
        if (!manufactured && maximum > 1e-7)
            throw std::runtime_error("Uniform state generated a nonzero RHS: " +
                                     std::to_string(maximum));
        if constexpr (std::is_same_v<Model, SingleDiffusiveField>)
            if (!manufactured)
                CheckParticleCoupling(this->m_session, this->m_graph,
                                      this->m_expansionLists[0]);
        if constexpr (std::is_same_v<Model, ReducedBraginskii>)
        {
            this->E.template Initialize<NektarSpaces::HostSpace>(0.);
            for (auto &block : this->E.GetBlocks())
            {
                auto electric =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                std::fill(electric + 2 * block.CompSize(),
                          electric + 3 * block.CompSize(), 1.);
            }
            this->m_math.zero(rhs);
            this->AddPressureWork(this->m_fields, rhs, factor);
            for (auto &block : rhs.GetBlocks())
            {
                auto force =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                const unsigned momentum = this->GetMomentumId(0);
                for (size_t e = 0; e < block.GetNumElements(); ++e)
                    for (size_t q = 0; q < block.GetNumData(); ++q)
                    {
                        const size_t width = block.GetInterleaveWidth();
                        const size_t p =
                            (e / width) * block.GetNumData() * width +
                            q * width + e % width;
                        if (std::abs(force[momentum * block.CompSize() + p] -
                                     2 * factor) > 1e-8)
                            throw std::runtime_error(
                                "Parallel electric-force sign or out-of-plane "
                                "component is incorrect.");
                    }
            }
            this->m_math.zero(this->E);
        }
        if constexpr (std::is_same_v<Model, ElectrostaticTurbulence>)
        {
            Field saved("Saved equilibrium", attrs, this->m_variables, 1);
            this->m_math.copy(this->m_fields, saved);
            coords.SetInterleaveWidth(this->m_fields);
            for (unsigned b = 0; b < this->m_fields.GetBlocks().size(); ++b)
            {
                auto &block = this->m_fields.GetBlocks()[b];
                auto u =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                auto xyz =
                    coords.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                const unsigned energy = this->GetEnergyId(0);
                for (size_t p = 0; p < block.CompSize(); ++p)
                    u[energy * block.CompSize() + p] += .1 * xyz[p];
            }
            this->PrepareAdvection(this->m_fields, 0.);
            this->m_extraVelocity.SetInterleaveWidth(this->m_fields);
            for (auto &block : this->m_extraVelocity.GetBlocks())
            {
                auto velocity =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                const size_t component = density * this->m_coordim + 1;
                for (size_t e = 0; e < block.GetNumElements(); ++e)
                    for (size_t q = 0; q < block.GetNumData(); ++q)
                    {
                        const size_t width = block.GetInterleaveWidth();
                        const size_t p =
                            (e / width) * block.GetNumData() * width +
                            q * width + e % width;
                        if (std::abs(
                                velocity[component * block.CompSize() + p]) >
                            1e-8)
                            throw std::runtime_error(
                                "ExB and ion diamagnetic drifts do not balance "
                                "in a pressure equilibrium.");
                    }
            }
            this->m_math.copy(saved, this->m_fields);
            this->PrepareAdvection(this->m_fields, 0.);
        }
        this->m_timeOp->Apply(this->m_fields);
        for (auto &block : this->m_fields.GetBlocks())
        {
            auto data =
                block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            for (size_t p = 0;
                 p < block.CompSize() * this->m_fields.GetNumComponents(); ++p)
                if (!std::isfinite(data[p]))
                    throw std::runtime_error(
                        "Explicit step produced a non-finite state.");
        }
        std::cout << "PASS: " << count << " points, maximum RHS " << maximum
                  << ", relative error "
                  << (reference > 0 ? std::sqrt(error / reference) : 0) << "\n";
    }
};
template <class Model>
void Run(const LibUtilities::SessionReaderSharedPtr &s,
         const SpatialDomains::MeshGraphSharedPtr &g, bool manufactured)
{
    std::make_shared<CheckSystem<Model>>(s, g)->Check(manufactured);
}
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int result = 0;
    try
    {
        CheckTransportCoefficients();
        auto session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        auto graph   = SpatialDomains::MeshGraphIO::Read(session);
        const auto model        = session->GetSolverInfo("EQTYPE");
        const bool manufactured = session->GetParameter("Manufactured") != 0;
        if (model == "SingleDiffusiveField")
            Run<SingleDiffusiveField>(session, graph, manufactured);
        else if (model == "DoubleDiffusiveField")
            Run<DoubleDiffusiveField>(session, graph, false);
        else if (model == "ReducedBraginskii")
            Run<ReducedBraginskii>(session, graph, false);
        else if (model == "ElectrostaticTurbulence")
            Run<ElectrostaticTurbulence>(session, graph, false);
        else if (model == "NullSystem")
            Run<NullSystem>(session, graph, false);
        else
            throw std::runtime_error("Unknown equation system.");
        session->Finalise();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << "\n";
        result = 1;
    }
    MPI_Finalize();
    return result;
}
