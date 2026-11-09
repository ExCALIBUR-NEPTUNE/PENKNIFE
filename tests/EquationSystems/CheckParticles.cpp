#include "ParticleSystems/ReactionSystem.hpp"
#include <MultiRegions/ElmtOps/Expression/ExpressionOp.hpp>
#include <stdexcept>
using namespace Nektar;
using namespace PENKNIFE;
class CouplingCheck : public ReactionSystem
{
public:
    CouplingCheck(NESOReaderSharedPtr config,
                  SpatialDomains::MeshGraphSharedPtr graph)
        : ReactionSystem(config, graph)
    {
    }
    void set_up_species() override
    {
    }
    void set_up_boundaries() override
    {
    }
    void pre_advection(ParticleSubGroupSharedPtr) override
    {
    }
    void Check(const MultiRegions::ExpListSharedPtr &expansion)
    {
        init_object();
        marker_group =
            std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);
        marker_cell_id_translation = std::make_shared<CellIDTranslation>(
            sycl_target, marker_group->cell_id_dat, particle_mesh_interface);
        auto add = [&](ParticleGroupSharedPtr group, double position)
        {
            ParticleSet particles(1, particle_spec);
            particles[Sym<REAL>("POSITION")][0][0]      = position;
            particles[Sym<REAL>("POSITION")][0][1]      = position;
            particles[Sym<INT>("CELL_ID")][0][0]        = 0;
            particles[Sym<INT>("INTERNAL_STATE")][0][0] = 0;
            group->add_particles_local(particles);
        };
        add(particle_group, -.5);
        add(marker_group, .5);
        transfer_particles();
        marker_group->hybrid_move();
        marker_cell_id_translation->execute();
        marker_group->cell_move();
        auto attrs = MultiRegions::GetBlockAttributes<double, FieldState::Phys>(
            expansion);
        PhysField values("Particle evaluation test", attrs, {"n"}, 1);
        values.Initialize<NektarSpaces::HostSpace>(0.);
        values.ReshapeStorage(2, "Serial");
        auto expression =
            MultiRegions::ExpressionOp<double>::Create(expansion, {"n"});
        expression->SetExpressions({std::make_shared<LibUtilities::Equation>(
            config->get_session()->GetInterpreter(), "2+x+y")});
        expression->Apply(values, values);
        setup_evaluate_fields(expansion, values,
                              {Sym<REAL>("ELECTRON_DENSITY")}, {0});
        evaluate_fields();
        auto copy = [&](ParticleGroupSharedPtr group)
        {
            particle_loop(
                group,
                [](auto density, auto weight) { weight.at(0) = density.at(0); },
                Access::read(Sym<REAL>("ELECTRON_DENSITY")),
                Access::write(Sym<REAL>("ELECTRON_SOURCE_ENERGY")))
                ->execute();
        };
        copy(particle_group);
        copy(marker_group);
        PhysField sources("Projection test", attrs, {"n"}, 1);
        CoeffField coefficients(
            "Projection test coefficients",
            MultiRegions::GetBlockAttributes<double, FieldState::Coeff>(
                expansion),
            {"n"}, 1);
        sources.Initialize<NektarSpaces::HostSpace>(0.);
        coefficients.Initialize<NektarSpaces::HostSpace>(0.);
        coefficients.ReshapeStorage(2, "Serial");
        const std::vector<Sym<REAL>> symbols{
            Sym<REAL>("ELECTRON_SOURCE_ENERGY")};
        const std::vector<int> components{0};
        ProjectTransformation ordinary(expansion, sources, coefficients,
                                       symbols, components, particle_group,
                                       cell_id_translation);
        ProjectTransformation markers(expansion, sources, coefficients, symbols,
                                      components, marker_group,
                                      marker_cell_id_translation);
        ordinary.transform(particle_sub_group(particle_group));
        ordinary.transform(particle_sub_group(particle_group));
        markers.transform(particle_sub_group(marker_group));
        sources.ReshapeStorage(1, "Serial");
        double integral  = 0;
        auto collections = MultiRegions::GetCollections(expansion);
        for (unsigned b = 0; b < collections.size(); ++b)
        {
            auto &block    = sources.GetBlocks()[b];
            auto data      = block.GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto &elements = collections[b].GetExpVector();
            for (unsigned e = 0; e < elements.size(); ++e)
            {
                Array<OneD, const double> value(block.GetNumData(),
                                                data + e * block.GetNumData());
                integral += elements[e]->Integral(value);
            }
        }
        if (std::abs(integral - 5.) > 1e-8)
            throw std::runtime_error("Particle evaluation/append/marker "
                                     "projection lost source weight: " +
                                     std::to_string(integral));
        free();
    }
};
void CheckParticleCoupling(const LibUtilities::SessionReaderSharedPtr &session,
                           const SpatialDomains::MeshGraphSharedPtr &graph,
                           const MultiRegions::ExpListSharedPtr &expansion)
{
    auto config = std::make_shared<NESOReader>(session);
    config->read_species();
    config->read_vantage();
    CouplingCheck particles(config, graph);
    particles.Check(expansion);
}
