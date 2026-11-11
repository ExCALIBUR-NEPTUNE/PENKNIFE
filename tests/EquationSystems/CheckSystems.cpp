#include "Closures/DiffusiveFields.hpp"
#include "EquationSystems/AnisotropicDiffusion.hpp"
#include "EquationSystems/ElectrostaticTurbulence.hpp"
#include "EquationSystems/NullSystem.hpp"
#include "EquationSystems/ReducedBraginskii.hpp"

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
    void CheckPointwiseParticles()
    {
        auto attrs = MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
            this->m_expansionLists[0]);
        const unsigned count = this->GetNumSpecies();
        auto op = PlasmaPointwiseOp<double>::Create(this->m_expansionLists[0],
                                                    this->m_variables);
        op->SetSpecies(this->n_ids, this->v_ids, this->e_ids, this->m_masses,
                       this->m_charges, this->GetIons(), this->ee_idx, this->me,
                       .25);
        std::vector<unsigned> density(count), momentum(count), energy(count);
        for (unsigned species = 0; species < count; ++species)
        {
            density[species]  = 5 * species;
            momentum[species] = 5 * species + 1;
            energy[species]   = 5 * species + 4;
        }
        op->SetSources(density, momentum, energy,
                       this->ee_idx == this->InvalidId ? this->InvalidId
                                                       : 5 * count);
        Field source("Synthetic particle sources", attrs, 5 * count + 1, 1);
        Field parallelCurrent("Test parallel current", attrs, 1, 1);
        Field quasineutralDensity("Test quasineutral density", attrs, 1, 1);
        Field quasineutralMomentum("Test quasineutral momentum", attrs, 1, 1);
        Field rhs("Particle RHS", attrs, this->m_variables, 1);
        Field electronDensity("Test electron density", attrs, 1, 1);
        Field electronMomentum("Test electron momentum", attrs, 1, 1);
        Field evaluated("Particle fields", attrs, 11 + 5 * count, 1);
        source.template Initialize<NektarSpaces::HostSpace>(3.);
        parallelCurrent.template Initialize<NektarSpaces::HostSpace>(3.);
        rhs.template Initialize<NektarSpaces::HostSpace>(.5);
        electronDensity.template Initialize<NektarSpaces::HostSpace>(2.);
        electronMomentum.template Initialize<NektarSpaces::HostSpace>(4. *
                                                                      this->me);
        op->AddSources(rhs, source, this->b_unit, .37, 1.25);
        op->ParticleFields(this->m_fields, electronDensity, electronMomentum,
                           this->b_unit, this->B, this->E, evaluated);
        op->Quasineutrality(this->m_fields, parallelCurrent,
                            quasineutralDensity, quasineutralMomentum);
        this->m_fields.SetInterleaveWidth(rhs);
        this->b_unit.SetInterleaveWidth(rhs);
        this->B.SetInterleaveWidth(rhs);
        this->E.SetInterleaveWidth(rhs);
        for (unsigned b = 0; b < rhs.GetBlocks().size(); ++b)
        {
            auto &block         = rhs.GetBlocks()[b];
            const size_t stride = block.CompSize(),
                         width  = block.GetInterleaveWidth(),
                         np     = block.GetNumData();
            auto actual =
                block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto values =
                evaluated.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto u = this->m_fields.GetBlocks()[b]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto unit =
                this->b_unit.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto magnetic =
                this->B.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto electric =
                this->E.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto current =
                parallelCurrent.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto quasineutralDensityValue =
                quasineutralDensity.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto quasineutralMomentumValue =
                quasineutralMomentum.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            for (size_t e = 0; e < block.GetNumElements(); ++e)
                for (size_t q = 0; q < np; ++q)
                {
                    const size_t p =
                        (e / width) * np * width + q * width + e % width;
                    auto check = [](double actual, double expected)
                    {
                        if (!std::isfinite(actual) ||
                            std::abs(actual - expected) >
                                1e-10 * (1 + std::abs(expected)))
                            throw std::runtime_error(
                                "Pointwise particle fields/source scaling is "
                                "incorrect.");
                    };
                    const double delta = .37 * 3. / 1.25;
                    for (unsigned c = 0; c < this->m_variables.size(); ++c)
                    {
                        double expected = .5;
                        if (c == this->ee_idx)
                            expected += delta;
                        for (unsigned species = 0; species < count; ++species)
                        {
                            if (c == this->n_ids[species] ||
                                c == this->e_ids[species])
                                expected += delta;
                            if (c == this->v_ids[species])
                                for (unsigned d = 0; d < 3; ++d)
                                    expected += delta * unit[d * stride + p];
                        }
                        check(actual[c * stride + p], expected);
                    }
                    check(values[p], 2.);
                    check(values[stride + p],
                          this->ee_idx == this->InvalidId
                              ? .25
                              : u[this->ee_idx * stride + p] / 3.);
                    for (unsigned d = 0; d < 3; ++d)
                    {
                        check(values[(2 + d) * stride + p],
                              2 * unit[d * stride + p]);
                        check(values[(5 + d) * stride + p],
                              electric[d * stride + p]);
                        check(values[(8 + d) * stride + p],
                              magnetic[d * stride + p]);
                    }
                    for (unsigned species = 0; species < count; ++species)
                    {
                        const unsigned c = 11 + 5 * species;
                        const double n =
                            this->HasDensity(species)
                                ? u[this->n_ids[species] * stride + p]
                                : 0.;
                        check(values[c * stride + p], n);
                        check(values[(c + 1) * stride + p],
                              n > 0 && this->HasEnergy(species)
                                  ? (2. / 3.) *
                                        u[this->e_ids[species] * stride + p] / n
                                  : .25);
                        for (unsigned d = 0; d < 3; ++d)
                            check(values[(c + 2 + d) * stride + p],
                                  n > 0 && this->HasMomentum(species)
                                      ? u[this->v_ids[species] * stride + p] *
                                            unit[d * stride + p] /
                                            (this->m_masses[species] * n)
                                      : 0.);
                    }
                    double expectedDensity  = 0.;
                    double expectedMomentum = -this->me * current[p];
                    for (const unsigned species : this->GetIons())
                    {
                        expectedDensity += this->m_charges[species] *
                                           u[this->n_ids[species] * stride + p];
                        if (this->HasMomentum(species))
                            expectedMomentum +=
                                (this->me / this->m_masses[species]) *
                                this->m_charges[species] *
                                u[this->v_ids[species] * stride + p];
                    }
                    check(quasineutralDensityValue[p], expectedDensity);
                    check(quasineutralMomentumValue[p], expectedMomentum);
                }
        }
    }
    void CheckPolarization(Field &coords)
    {
        auto &velocity = this->m_velocity;
        coords.SetInterleaveWidth(velocity);
        this->m_fields.SetInterleaveWidth(velocity);
        this->B.SetInterleaveWidth(velocity);
        this->mag_B.SetInterleaveWidth(velocity);
        for (unsigned b = 0; b < velocity.GetBlocks().size(); ++b)
        {
            auto &block = velocity.GetBlocks()[b];
            auto v =
                block.template GetPtr<NektarSpaces::HostSpace, ReadWrite>();
            auto x = coords.GetBlocks()[b]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            const size_t stride = block.CompSize();
            for (unsigned species : this->GetIons())
                for (unsigned d = 0; d < this->m_coordim; ++d)
                    for (size_t p = 0; p < stride; ++p)
                        v[(species * this->m_coordim + d) * stride + p] =
                            (species + 1.) * (d + 1.) *
                            (1. + .2 * std::sin(M_PI * x[d * stride + p]));
        }
        this->m_pointwise->PolarizationTensor(this->m_fields, velocity, this->B,
                                              this->mag_B, this->m_tensor,
                                              this->omega_c);
        this->m_tensorDerivative->Apply(this->m_tensor, this->m_tensorGradient);
        this->m_pointwise->PolarizationFlux(
            this->m_fields, this->m_tensorGradient, this->m_omegaFlux);
        Field reference(
            "Reference polarization flux",
            MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
                this->m_expansionLists[0]),
            this->m_coordim, 1);
        reference.SetInterleaveWidth(this->m_omegaFlux);
        this->m_math.zero(reference);
        for (unsigned species : this->GetIons())
        {
            this->m_math.zero(this->m_tensor);
            this->m_tensor.SetInterleaveWidth(velocity);
            for (unsigned b = 0; b < velocity.GetBlocks().size(); ++b)
            {
                auto &block         = velocity.GetBlocks()[b];
                const size_t stride = block.CompSize();
                auto v =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                auto u =
                    this->m_fields.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                auto B =
                    this->B.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                auto B2 =
                    this->mag_B.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                auto tensor =
                    this->m_tensor.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                for (size_t p = 0; p < stride; ++p)
                {
                    if (B2[p] <= 0)
                        continue;
                    double pointVelocity[3] = {};
                    for (unsigned d = 0; d < this->m_coordim; ++d)
                        pointVelocity[d] =
                            v[(species * this->m_coordim + d) * stride + p];
                    for (unsigned d = 0; d < this->m_coordim; ++d)
                    {
                        const unsigned j = (d + 1) % 3, k = (d + 2) % 3;
                        const double w =
                            u[this->n_ids[species] * stride + p] *
                            (pointVelocity[j] * B[k * stride + p] -
                             pointVelocity[k] * B[j * stride + p]) /
                            B2[p];
                        for (unsigned t = 0; t < this->m_coordim; ++t)
                            tensor[(t * this->m_coordim + d) * stride + p] =
                                w * pointVelocity[t];
                    }
                }
            }
            this->m_tensorDerivative->Apply(this->m_tensor,
                                            this->m_tensorGradient);
            this->m_tensorGradient.SetInterleaveWidth(reference);
            for (unsigned b = 0; b < reference.GetBlocks().size(); ++b)
            {
                auto &block         = reference.GetBlocks()[b];
                const size_t stride = block.CompSize();
                auto expected =
                    block.template GetPtr<NektarSpaces::HostSpace, ReadWrite>();
                auto derivative =
                    this->m_tensorGradient.GetBlocks()[b]
                        .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
                for (size_t p = 0; p < stride; ++p)
                    for (unsigned d = 0; d < this->m_coordim; ++d)
                        for (unsigned j = 0; j < this->m_coordim; ++j)
                            expected[d * stride + p] +=
                                this->m_masses[species] / this->omega_c *
                                derivative[((d * this->m_coordim + j) *
                                                this->m_coordim +
                                            j) *
                                               stride +
                                           p];
            }
        }
        this->m_omegaFlux.SetInterleaveWidth(reference);
        double maximum = 0;
        for (unsigned b = 0; b < reference.GetBlocks().size(); ++b)
        {
            auto &block         = reference.GetBlocks()[b];
            const size_t stride = block.CompSize(),
                         width  = block.GetInterleaveWidth(),
                         np     = block.GetNumData();
            auto expected =
                block.template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto actual =
                this->m_omegaFlux.GetBlocks()[b]
                    .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            for (size_t e = 0; e < block.GetNumElements(); ++e)
                for (size_t q = 0; q < np; ++q)
                    for (unsigned d = 0; d < this->m_coordim; ++d)
                    {
                        const size_t p = d * stride + (e / width) * np * width +
                                         q * width + e % width;
                        maximum = std::max(maximum, std::abs(expected[p]));
                        if (std::abs(actual[p] - expected[p]) >
                            1e-9 * (1 + std::abs(expected[p])))
                            throw std::runtime_error(
                                "Batched polarization differs from per-species "
                                "derivatives.");
                    }
        }
        if (maximum < 1e-8)
            throw std::runtime_error(
                "Polarization regression has a zero reference flux.");
    }
    void Check(bool manufactured)
    {
        this->InitObject();
        this->SetInitialConditions();
        CheckPointwiseParticles();
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
        if constexpr (std::is_same_v<Model, AnisotropicDiffusion>)
            if (!manufactured && this->DiffusionModel() == "SingleDiffusive")
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
                for (unsigned species : this->GetIons())
                {
                    const unsigned energy = this->GetEnergyId(species);
                    for (size_t p = 0; p < block.CompSize(); ++p)
                        u[energy * block.CompSize() + p] +=
                            .1 * this->m_charges[species] * xyz[p];
                }
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
            CheckPolarization(coords);
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
        if (model == "AnisotropicDiffusion")
            Run<AnisotropicDiffusion>(session, graph, manufactured);
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
