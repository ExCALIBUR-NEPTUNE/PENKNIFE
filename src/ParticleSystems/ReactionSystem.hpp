#ifndef REACTIONSYSTEM_HPP
#define REACTIONSYSTEM_HPP

#include <reactions/reactions.hpp>

#include "ParticleSystem.hpp"

using namespace VANTAGE::Reactions;

namespace PENKNIFE
{

class ReactionSystem : public ParticleSystem
{

public:
    static std::string class_name;
    static std::shared_ptr<ParticleSystem> create(
        const NESOReaderSharedPtr session, const SD::MeshGraphSharedPtr graph)
    {
        std::shared_ptr<ParticleSystem> p =
            MemoryManager<ReactionSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    ReactionSystem(NESOReaderSharedPtr session, SD::MeshGraphSharedPtr graph);

    ~ReactionSystem() override = default;

    void free() override;

    inline void evaluate_fields() override
    {
        ParticleSystem::evaluate_fields();
        if (this->marker_group)
        {
            if (!marker_evaluate)
            {
                marker_evaluate = FieldEvaluateOp<double>::Create(
                    proto_field, eval_fields->GetComponentNames());
                marker_evaluate->SetParticles(marker_group,
                                              marker_cell_id_translation);
            }
            for (unsigned c = 0; c < eval_syms.size(); ++c)
                marker_evaluate->Apply(*eval_fields, eval_syms[c], c,
                                       eval_comps[c]);
            set_marker_weights();
        }
    }

    inline void apply_timestep(const double dt) override
    {
        this->zeroer_transform->transform(
            particle_sub_group(this->particle_group));
        if (this->marker_group)
            this->zeroer_transform->transform(
                particle_sub_group(this->marker_group));

        ParticleSystem::apply_timestep(dt);

        if (reaction_controller)
            reaction_controller->apply(this->particle_group, dt);
        if (recomb_controller)
            recomb_controller->apply(this->marker_group, dt,
                                     this->particle_group);

        auto partitions = particle_group_partition(this->particle_group,
                                                   Sym<INT>("INTERNAL_STATE"),
                                                   this->species_map.size());

        int s = 0;
        for (const auto &[k, v] : this->species_map)
        {
            species_map[k].sub_group = partitions[s++];
        };
    }

    void set_up_reactions();
    void set_up_boundaries() override;

    inline void pre_advection(ParticleSubGroupSharedPtr sg) override
    {
        this->boundary->pre_advection(sg);
    };

    inline void apply_boundary_conditions(ParticleSubGroupSharedPtr sg,
                                          ParticleGroupSharedPtr cg,
                                          double dt) override
    {
        this->boundary->execute(sg, cg, dt);
    };

    void finish_setup(const std::shared_ptr<ExpList> &expansion,
                      PhysField &src_fields, CoeffField &src_coeffs,
                      std::vector<Sym<REAL>> &syms,
                      std::vector<int> &components) override;

    void setup_reaction_controller();
    void setup_recomb_controller();

    void diag_setup() override;
    void diag_project() override;

    void append_diagnostics(
        Nektar::SolverCore::CoefficientFieldViews &views) override
    {
        ParticleSystem::append_diagnostics(views);
        for (auto &[name, species] : marker_map)
            views.push_back({marker_diag_coeffs.at(species.id),
                             proto_field,
                             unsigned(ndim),
                             {{name + "_MARKER_DENSITY", 0}}});
    }

    void create_markers(std::string k);
    void set_marker_weights();

    void output_setup(std::vector<Sym<REAL>> &syms) override;

    void write(const int step) override
    {
        ParticleSystem::write(step);
        if (this->h5part_marker)
        {
            this->h5part_marker->write();
        }
    }

    class ReactionsBoundary
    {

    public:
        ReactionsBoundary(
            Sym<REAL> time_step_prop_sym, SYCLTargetSharedPtr sycl_target,
            std::shared_ptr<ParticleMeshInterface> mesh,
            NESOReaderSharedPtr config,
            std::map<std::string, SpeciesInfo> &species,
            ParameterStoreSharedPtr store = std::make_shared<ParameterStore>());

        inline void pre_advection(ParticleSubGroupSharedPtr particle_sub_group)
        {
            this->composite_intersection->pre_integration(particle_sub_group);
        }

        inline void execute(ParticleSubGroupSharedPtr particle_sub_group,
                            ParticleGroupSharedPtr child_group, double dt)
        {
            NESOASSERT(this->ndim == 3 || this->ndim == 2,
                       "Unexpected number of dimensions.");
            auto groups = this->composite_intersection->get_intersections(
                particle_sub_group);

            for (auto &[id, sg] : groups)
            {
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<REAL>("NESO_PARTICLES_BOUNDARY_INTERSECTION_POINT"),
                    Sym<REAL>("NESO_PARTICLES_BOUNDARY_INTERSECTION_POINT"));
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<REAL>("NESO_PARTICLES_BOUNDARY_NORMAL"),
                    Sym<REAL>("NESO_PARTICLES_BOUNDARY_NORMAL"));
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<INT>("NESO_PARTICLES_BOUNDARY_METADATA"),
                    Sym<INT>("NESO_PARTICLES_BOUNDARY_METADATA"));

                this->boundary_truncation->execute(
                    sg,
                    get_particle_group(particle_sub_group)->position_dat->sym,
                    this->time_step_prop_sym,
                    this->composite_intersection->previous_position_sym);
                this->reaction_controllers[id]->apply(
                    sg, dt, child_group, ControllerMode::surface_mode);
            }
            remove_wrapper->transform(particle_sub_group);
        }

    private:
        Sym<REAL> time_step_prop_sym;

        SYCLTargetSharedPtr sycl_target;
        std::shared_ptr<CompositeInteraction::CompositeIntersection>
            composite_intersection;
        std::shared_ptr<BoundaryTruncation> boundary_truncation;

        std::map<int, std::shared_ptr<ReactionController>> reaction_controllers;
        std::shared_ptr<TransformationWrapper> remove_wrapper;

        const int ndim;
        const int vdim;
        REAL reset_distance;

        NESOReaderSharedPtr config;
    };

protected:
    class ProjectTransformation : public TransformationStrategy
    {

    public:
        ProjectTransformation(const std::shared_ptr<ExpList> &expansion,
                              PhysField &fields, CoeffField &coeffs,
                              const std::vector<Sym<REAL>> &syms,
                              const std::vector<int> &components,
                              ParticleGroupSharedPtr group,
                              CellIDTranslationSharedPtr translation)
            : fields(fields), coeffs(coeffs), syms(syms), components(components)
        {
            field_project = FieldProjectOp<double>::Create(
                expansion, fields.GetComponentNames());
            field_project->SetParticles(group, translation);
            field_project->SetAppend(true);
        }

        void transform(ParticleSubGroupSharedPtr group) override
        {
            field_project->Apply(group, syms, components, coeffs, fields);
        }

    private:
        PhysField &fields;
        CoeffField &coeffs;
        std::vector<Sym<REAL>> syms;
        std::vector<int> components;
        std::shared_ptr<FieldProjectOp<double>> field_project;
    };
    PhysField *projection_fields  = nullptr;
    CoeffField *projection_coeffs = nullptr;
    std::shared_ptr<TransformationWrapper> marker_project_wrapper;
    std::shared_ptr<FieldEvaluateOp<double>> marker_evaluate;

    std::shared_ptr<TransformationStrategy> project_transform;
    std::shared_ptr<TransformationStrategy> remove_transform;
    std::shared_ptr<TransformationStrategy> merge_transform;
    std::shared_ptr<TransformationWrapper> remove_wrapper;
    std::shared_ptr<TransformationWrapper> project_wrapper;
    std::shared_ptr<TransformationWrapper> merge_wrapper;
    std::shared_ptr<TransformationStrategy> zeroer_transform;
    std::shared_ptr<TransformationWrapper> zeroer_wrapper;

    uint64_t total_num_markers_added = 0;
    ParticleGroupSharedPtr marker_group;
    std::shared_ptr<CellIDTranslation> marker_cell_id_translation;

    std::map<std::string, SpeciesInfo> marker_map;
    std::shared_ptr<HostPerParticleBlockRNG<REAL>> rng_kernel;
    /// Reaction Controller
    std::shared_ptr<ReactionController> reaction_controller;
    std::shared_ptr<ReactionController> recomb_controller;

    std::vector<Sym<REAL>> marker_diag_syms;
    std::vector<int> marker_diag_components;
    std::shared_ptr<FieldProjectOp<double>> marker_diagnostic_project;
    std::map<int, PhysField> marker_diag_fields;
    std::map<int, CoeffField> marker_diag_coeffs;
    std::shared_ptr<H5Part> h5part_marker;

    std::shared_ptr<ReactionsBoundary> boundary;
};

} // namespace PENKNIFE
#endif