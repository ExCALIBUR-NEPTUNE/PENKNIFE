#ifndef REACTIONSYSTEM_HPP
#define REACTIONSYSTEM_HPP

#include "AMJUEL.hpp"
#include "ParticleSystem.hpp"
#include <neso_rng_toolkit.hpp>

using namespace VANTAGE::Reactions;

namespace NESO::Solvers::tokamak
{

class ReactionSystem : public ParticleSystem
{

public:
    static std::string class_name;
    static ParticleSystemSharedPtr create(const NESOReaderSharedPtr session,
                                          const SD::MeshGraphSharedPtr graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ReactionSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    ReactionSystem(NESOReaderSharedPtr session, SD::MeshGraphSharedPtr graph);

    ~ReactionSystem() override = default;

    std::shared_ptr<TransformationWrapper> zeroer_transform_wrapper;

    inline void integrate(const double time_end, const double dt) override
    {
        this->zeroer_transform_wrapper->transform(this->particle_group);
        ParticleSystem::integrate(time_end, dt);
        this->field_project->project(this->particle_group, this->src_syms,
                                     this->src_components);
    }

    inline void set_up_reactions()
    {
        auto prop_map = get_default_map();

        auto sycl_target = particle_group->sycl_target;
        const int rank   = sycl_target->comm_pair.rank_parent;

        std::uint64_t root_seed = 141351;
        std::uint64_t seed      = NESO::RNGToolkit::create_seeds(
            sycl_target->comm_pair.size_parent, rank, root_seed);

        // Create a Uniform distribution
        auto rng_normal = NESO::RNGToolkit::create_rng<REAL>(
            NESO::RNGToolkit::Distribution::Uniform<REAL>{
                NESO::RNGToolkit::Distribution::next_value(0.0), 1.0},
            seed, sycl_target->device, sycl_target->device_index);

        // Create an interface between NESO-RNG-Toolkit and NESO-Particles
        // KernelRNG
        auto rng_interface =
            make_rng_generation_function<GenericDeviceRNGGenerationFunction,
                                         REAL>(
                [=](REAL *d_ptr, const std::size_t num_samples) -> int
                { return rng_normal->get_samples(d_ptr, num_samples); });

        auto rng_kernel =
            host_atomic_block_kernel_rng<REAL>(rng_interface, 4 * 10);

        for (const auto &v : this->config->get_reactions())
        {
            std::shared_ptr<AbstractReaction> reaction;
            if (std::get<0>(v) == "Ionisation")
            {
                auto target_species =
                    Species(std::get<1>(v)[0],
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            this->species_map[std::get<1>(v)[0]].id);

                auto electron_species = Species("ELECTRON");

                if (std::get<2>(v).first == "Fixed")
                {
                    auto ionise_data        = FixedRateData(std::get<2>(v).second);
                    auto ionise_energy_data = FixedRateData(std::get<2>(v).second);
                    if (this->ndim == 2)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            2>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                    else if (this->ndim == 3)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            3>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto ionise_rate_data   = AMJUEL::ionise_rate_data();
                    auto ionise_energy_data = AMJUEL::ionise_energy_data();

                    if (this->ndim == 2)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_rate_data),
                            decltype(ionise_energy_data), 2>>(
                            this->particle_group->sycl_target, ionise_rate_data,
                            ionise_energy_data, target_species,
                            electron_species);
                    }
                    else if (this->ndim == 3)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_rate_data),
                            decltype(ionise_energy_data), 3>>(
                            this->particle_group->sycl_target, ionise_rate_data,
                            ionise_energy_data, target_species,
                            electron_species);
                    }
                }
            }
            else if (std::get<0>(v) == "Recombination")
            {
                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto neutral_species =
                    Species(std::get<1>(v)[0],
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            this->species_map[std::get<1>(v)[0]].id);

                auto marker_species =
                    Species(std::get<1>(v)[0],
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge - 1,
                            -1 - this->species_map[std::get<1>(v)[0]].id);

                if (std::get<2>(v).first == "Fixed")
                {
                    auto recomb_data = FixedRateData(1.0);
                    auto data1       = FixedRateData(1.0);
                    auto data2       = FixedRateData(1.0);
                    auto data_calculator =
                        DataCalculator<FixedRateData, FixedRateData,
                                       FixedRateData>(data1, data1, data2);
                    if (this->ndim == 2)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<2>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<3>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto recomb_rate_data   = AMJUEL::recomb_rate_data();
                    auto recomb_energy_data = AMJUEL::recomb_energy_data();

                    auto constant_rate_cross_section =
                        ConstantRateCrossSection(1.0);
                    auto recomb_data_calc_sampler = FilteredMaxwellianSampler<
                        2, decltype(constant_rate_cross_section)>(
                        (constants::temp_SI * constants::k_B) /
                            (marker_species.get_mass() * norm::mass_amu_SI *
                             norm::vel * norm::vel),
                        constant_rate_cross_section, rng_kernel);
                    auto data_calculator =
                        DataCalculator<decltype(recomb_energy_data),
                                       decltype(recomb_data_calc_sampler)>(
                            recomb_energy_data, recomb_data_calc_sampler);

                    if (this->ndim == 2)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<2>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_rate_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_rate_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<3>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_rate_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_rate_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                }
            }

            else if (std::get<0>(v) == "ChargeExchange")
            {
                auto projectile_species =
                    Species(std::get<1>(v)[0],
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            this->species_map[std::get<1>(v)[0]].id);
                auto target_species =
                    Species(std::get<1>(v)[1],
                            this->species_map[std::get<1>(v)[1]].mass,
                            this->species_map[std::get<1>(v)[1]].charge,
                            this->species_map[std::get<1>(v)[1]].id);
                auto parent_mass = target_species.get_mass();
                auto child_mass  = projectile_species.get_mass();
                auto reduced_mass =
                    (parent_mass * child_mass) / (parent_mass + child_mass);

                if (std::get<2>(v).first == "Fixed")
                {
                    auto rate_data              = FixedRateData(1.0);
                    auto constant_cross_section = ConstantRateCrossSection(1.0);

                    auto data_calc_sampler = FilteredMaxwellianSampler<
                        2, decltype(constant_cross_section)>(
                        (constants::temp_SI * constants::k_B) /
                            (child_mass * norm::mass_amu_SI * norm::vel *
                             norm::vel),
                        constant_cross_section, rng_kernel);
                    auto data_calculator =
                        DataCalculator<decltype(data_calc_sampler)>(
                            data_calc_sampler);
                    if (this->ndim == 2)
                    {
                        auto cx_kernel = CXReactionKernels<2>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto cx_kernel = CXReactionKernels<3>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto rate_data =
                        AMJUEL::cx_rate_data(parent_mass, child_mass);
                    auto cross_section =
                        AMJUEL::amjuel_fit_cross_section(reduced_mass);

                    auto data_calc_sampler =
                        FilteredMaxwellianSampler<2, decltype(cross_section)>(
                            (constants::temp_SI * constants::k_B) /
                                (child_mass * norm::mass_amu_SI * norm::vel *
                                 norm::vel),
                            cross_section, rng_kernel);

                    auto data_calculator =
                        DataCalculator<decltype(data_calc_sampler)>(
                            data_calc_sampler);

                    if (this->ndim == 2)
                    {
                        auto cx_kernel = CXReactionKernels<2>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto cx_kernel = CXReactionKernels<3>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                }
            }
            this->reaction_controller->add_reaction(reaction);
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
            ParameterStoreSharedPtr store = std::make_shared<ParameterStore>())
            : time_step_prop_sym(time_step_prop_sym), sycl_target(sycl_target),
              ndim(mesh->get_ndim()), config(config)
        {
            auto reflection_removal_wrapper =
                std::make_shared<TransformationWrapper>(
                    std::vector<std::shared_ptr<MarkingStrategy>>{
                        make_marking_strategy<
                            ComparisonMarkerSingle<REAL, LessThanComp>>(
                            Sym<REAL>("WEIGHT"), 1.0e-12)},
                    make_transformation_strategy<
                        SimpleRemovalTransformationStrategy>());
            config->read_boundary_regions();

            for (auto &v : this->config->get_surface_reactions())
            {
                std::vector<int> boundary_ids = std::get<2>(v);
                for (int b_id : boundary_ids)
                {
                    if (!this->reaction_controllers[b_id])
                    {
                        reaction_controllers[b_id] =
                            std::make_shared<ReactionController>(
                                std::vector<
                                    std::shared_ptr<TransformationWrapper>>{},
                                std::vector<
                                    std::shared_ptr<TransformationWrapper>>{});
                    }

                    if (std::get<0>(v) == "Specular")
                    {
                        for (const auto &s : std::get<1>(v))
                        {
                            int s_id       = species[s].id;
                            auto rate_data = FixedRateData(1.0);
                            if (this->ndim == 2)
                            {
                                auto reflection_kernels =
                                    SpecularReflectionKernels<2>();

                                auto reflection_reaction =
                                    std::make_shared<LinearReactionBase<
                                        0, FixedRateData,
                                        SpecularReflectionKernels<2>>>(
                                        sycl_target, s_id, std::array<int, 0>{},
                                        rate_data, reflection_kernels);

                                this->reaction_controllers[b_id]->add_reaction(
                                    reflection_reaction);
                            }
                            else if (this->ndim == 3)
                            {
                                auto reflection_kernels =
                                    SpecularReflectionKernels<3>();

                                auto reflection_reaction =
                                    std::make_shared<LinearReactionBase<
                                        0, FixedRateData,
                                        SpecularReflectionKernels<3>>>(
                                        sycl_target, s_id, std::array<int, 0>{},
                                        rate_data, reflection_kernels);

                                this->reaction_controllers[b_id]->add_reaction(
                                    reflection_reaction);
                            }
                        }
                    }
                }
            }

            this->composite_intersection =
                std::make_shared<CompositeInteraction::CompositeIntersection>(
                    this->sycl_target, mesh, config->get_boundary_regions());

            this->reset_distance =
                store->get<REAL>("ReactionsBoundary/reset_distance", 1.0e-4);

            this->boundary_truncation = std::make_shared<BoundaryTruncation>(
                this->ndim, this->reset_distance);
        }

        void pre_advection(ParticleSubGroupSharedPtr particle_sub_group)
        {
            this->composite_intersection->pre_integration(particle_sub_group);
        }

        void execute(ParticleSubGroupSharedPtr particle_sub_group, double dt)
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
                    sg, dt, ControllerMode::surface_mode);
            }
        }

    private:
        Sym<REAL> time_step_prop_sym;

        SYCLTargetSharedPtr sycl_target;
        std::shared_ptr<CompositeInteraction::CompositeIntersection>
            composite_intersection;
        std::shared_ptr<BoundaryTruncation> boundary_truncation;

        std::map<int, std::shared_ptr<ReactionController>> reaction_controllers;

        int ndim;
        REAL reset_distance;

        NESOReaderSharedPtr config;
    };

    void set_up_boundaries() override;

    void pre_advection(ParticleSubGroupSharedPtr sg) override
    {
        this->boundary->pre_advection(sg);
    };

    void apply_boundary_conditions(ParticleSubGroupSharedPtr sg,
                                   double dt) override
    {
        this->boundary->execute(sg, dt);
    };

    inline void integrate_inner_ion(ParticleSubGroupSharedPtr sg,
                                    const double dt_inner) override
    {
        ParticleSystem::integrate_inner_ion(sg, dt_inner);
        if (this->config->get_reactions().size())
            reaction_controller->apply(sg, dt_inner);
    }
    inline void integrate_inner_neutral(ParticleSubGroupSharedPtr sg,
                                        const double dt_inner) override
    {
        ParticleSystem::integrate_inner_neutral(sg, dt_inner);
        if (this->config->get_reactions().size())
            reaction_controller->apply(sg, dt_inner);
    }

    inline void finish_setup(
        std::vector<std::shared_ptr<DisContField>> &src_fields,
        std::vector<Sym<REAL>> &syms, std::vector<int> &components) override
    {
        this->src_syms       = syms;
        this->src_components = components;

        this->field_project = std::make_shared<FieldProject<DisContField>>(
            src_fields, this->particle_group, this->cell_id_translation);

        auto project_transform = std::make_shared<ProjectTransformation>(
            src_fields, this->src_syms, this->src_components,
            this->particle_group, this->cell_id_translation);
        auto project_transform_wrapper =
            std::make_shared<TransformationWrapper>(
                std::dynamic_pointer_cast<TransformationStrategy>(
                    project_transform));

        auto remove_transform =
            std::make_shared<SimpleRemovalTransformationStrategy>();
        auto remove_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            1e-10)},
            std::dynamic_pointer_cast<TransformationStrategy>(
                remove_transform));

        std::shared_ptr<TransformationStrategy> merge_transform;

        if (this->ndim == 2)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<2>>();
        }
        else if (this->ndim == 3)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<3>>();
        }

        auto merge_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            0.0000000001)},
            merge_transform);

        std::vector<std::string> src_names{"ELECTRON_SOURCE_DENSITY",
                                           "ELECTRON_SOURCE_ENERGY",
                                           "ELECTRON_SOURCE_MOMENTUM"};

        for (auto &[k, v] : this->config->get_particle_species())
        {
            src_names.push_back(k + "_SOURCE_DENSITY");
            src_names.push_back(k + "_SOURCE_ENERGY");
            src_names.push_back(k + "_SOURCE_MOMENTUM");
        }

        auto zeroer_transform =
            make_transformation_strategy<ParticleDatZeroer<REAL>>(src_names);

        this->zeroer_transform_wrapper =
            std::make_shared<TransformationWrapper>(
                std::dynamic_pointer_cast<TransformationStrategy>(
                    zeroer_transform));

        this->reaction_controller = std::make_shared<ReactionController>(
            std::vector<std::shared_ptr<TransformationWrapper>>{
                /*zeroer_transform_wrapper, merge_transform_wrapper,
                remove_transform_wrapper*/},
            std::vector<std::shared_ptr<TransformationWrapper>>{
                /*project_transform_wrapper, merge_transform_wrapper,
                remove_transform_wrapper*/});

        set_up_reactions();
        init_output("particle_trajectory.h5part", Sym<REAL>("POSITION"),
                    Sym<INT>("INTERNAL_STATE"), Sym<INT>("CELL_ID"),
                    Sym<REAL>("VELOCITY"), Sym<REAL>("MAGNETIC_FIELD"),
                    Sym<REAL>("ELECTRON_DENSITY"), this->src_syms,
                    Sym<REAL>("WEIGHT"), Sym<INT>("ID"),
                    Sym<REAL>("TOT_REACTION_RATE"), Sym<REAL>("FLUID_DENSITY"),
                    Sym<REAL>("FLUID_TEMPERATURE"));
    }

protected:
    class ProjectTransformation : public TransformationStrategy
    {

    public:
        ProjectTransformation(
            std::vector<std::shared_ptr<DisContField>> &src_fields,
            std::vector<Sym<REAL>> &src_syms, std::vector<int> &src_components,
            ParticleGroupSharedPtr particle_group,
            std::shared_ptr<CellIDTranslation> cell_id_translation)
            : syms(src_syms), components(src_components)
        {
            this->field_project = std::make_shared<FieldProject<DisContField>>(
                src_fields, particle_group, cell_id_translation);
        }

        void transform(ParticleSubGroupSharedPtr sub_group) override
        {
            this->field_project->project(sub_group, syms, components);
        }

    private:
        std::vector<Sym<REAL>> syms;
        std::vector<int> components;

        std::shared_ptr<FieldProject<DisContField>> field_project;
    };

    /// Reaction Controller
    std::shared_ptr<ReactionController> reaction_controller;

    std::shared_ptr<ReactionsBoundary> boundary;
};

} // namespace NESO::Solvers::tokamak
#endif