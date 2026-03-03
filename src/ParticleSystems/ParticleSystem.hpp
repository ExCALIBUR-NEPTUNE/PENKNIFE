#ifndef PENKNIFE_PARTICLE_SYSTEM_H
#define PENKNIFE_PARTICLE_SYSTEM_H

#include <array>

#include "../Misc/Constants.hpp"
#include <nektar_interface/function_evaluation.hpp>
#include <nektar_interface/function_projection.hpp>
#include <nektar_interface/particle_boundary_conditions.hpp>
#include <nektar_interface/particle_cell_mapping/particle_cell_mapping_common.hpp>
#include <nektar_interface/solver_base/partsys_base.hpp>
#include <nektar_interface/utilities.hpp>
#include <neso_particles.hpp>

namespace PENKNIFE
{

/**
 * @brief
 */
class ParticleSystem : public PartSysBase
{

public:
    static std::string class_name;
    /**
     * @brief Create an instance of this class and initialise it.
     */
    static ParticleSystemSharedPtr create(const NESOReaderSharedPtr &session,
                                          const SD::MeshGraphSharedPtr &graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ParticleSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    /**
     *  Create a new instance.
     *
     *  @param session Particle reader to use for parameters and simulation
     * specification.
     *  @param graph Nektar++ MeshGraph on which particles exist.
     *  @param comm (optional) MPI communicator to use - default MPI_COMM_WORLD.
     *
     */
    ParticleSystem(NESOReaderSharedPtr session, SD::MeshGraphSharedPtr graph,
                   MPI_Comm comm = MPI_COMM_WORLD);

    virtual ~ParticleSystem() override;

    /// Disable (implicit) copies.
    ParticleSystem(const ParticleSystem &st) = delete;
    /// Disable (implicit) copies.
    ParticleSystem &operator=(ParticleSystem const &a) = delete;

    virtual void init_spec() override;
    virtual void init_object() override;
    virtual void set_up_species() override;
    virtual void set_up_boundaries();

    struct SpeciesInfo
    {
        int id;
        double mass;
        double charge;

        std::shared_ptr<ParticleSubGroup> sub_group;
    };
    virtual std::map<std::string, SpeciesInfo> &get_species()
    {
        return species_map;
    }

    /**
     *  Integrate the particle system forward to the requested time using
     *  (at most) the requested time step.
     *
     *  @param time_end Target time to integrate to.
     *  @param dt Time step size.
     */
    inline virtual void integrate(const double time_end, const double dt)
    {
        // Get the current simulation time.
        NESOASSERT(time_end >= this->simulation_time,
                   "Cannot integrate backwards in time.");
        if (time_end == this->simulation_time)
        {
            return;
        }

        double time_tmp = this->simulation_time;
        while (time_tmp < time_end)
        {
            const double dt_inner = std::min(dt, time_end - time_tmp);
            this->add_sources(time_tmp, dt_inner);
            this->add_sinks(time_tmp, dt_inner);
            this->apply_timestep(dt_inner);
            this->transfer_particles();

            time_tmp += dt_inner;
        }

        this->simulation_time = time_end;
    }

    /**
     * Setup the projection object to use the src fields from the eqnsys.
     *
     * @param src_fields Nektar++ fields to project ionised particle data onto.
     * @param syms Corresponding Particle Syms
     * @param syms Corresponding components
     */
    virtual void finish_setup(
        std::vector<std::shared_ptr<DisContField>> &src_fields,
        std::vector<Sym<REAL>> &syms, std::vector<int> &components);

    virtual void diag_setup(const std::shared_ptr<DisContField> &diag_field);

    inline virtual void diag_project()
    {
        std::vector<Sym<REAL>> syms{Sym<REAL>("WEIGHT")};
        std::vector<int> components{0};
        this->diagnostic_project->project(this->particle_group, syms,
                                          components);
    }

    void add_sources(double time, double dt);
    void add_sinks(double time, double dt);

    /**
     *  Project the plasma source and momentum contributions from particle data
     *  onto field data.
     */
    inline virtual void project_source_terms()
    {
        NESOASSERT(this->field_project != nullptr,
                   "Field project object is null. Was setup_project called?");

        this->field_project->project(this->particle_group, this->src_syms,
                                     this->src_components);

        // remove fully ionised particles from the simulation
        remove_marked_particles();
    }

    virtual void setup_evaluate_fields(
        Array<OneD, std::shared_ptr<DisContField>> &E,
        Array<OneD, std::shared_ptr<DisContField>> &B,
        std::shared_ptr<DisContField> ne, std::shared_ptr<DisContField> Te,
        Array<OneD, std::shared_ptr<DisContField>> &ve);

    /**
     * Evaluate E and B at the particle locations.
     */
    inline virtual void evaluate_fields(
        Array<OneD, std::shared_ptr<DisContField>> &E,
        Array<OneD, std::shared_ptr<DisContField>> &B,
        std::shared_ptr<DisContField> ne, std::shared_ptr<DisContField> Te,
        Array<OneD, std::shared_ptr<DisContField>> &ve)
    {

        for (int d = 0; d < this->ndim; ++d)
        {
            NESOASSERT(this->field_evaluate_E[d] != nullptr,
                       "FieldEvaluate not setup.");
            this->field_evaluate_E[d]->evaluate(this->particle_group,
                                                Sym<REAL>("ELECTRIC_FIELD"), d,
                                                E[d]->GetCoeffs());
            NESOASSERT(this->field_evaluate_B[d] != nullptr,
                       "FieldEvaluate not setup.");
            this->field_evaluate_B[d]->evaluate(this->particle_group,
                                                Sym<REAL>("MAGNETIC_FIELD"), d,
                                                B[d]->GetCoeffs());
        }

        NESOASSERT(this->field_evaluate_ne != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_ne->evaluate(this->particle_group,
                                          Sym<REAL>("FLUID_DENSITY"), 0,
                                          ne->GetCoeffs());
        if (field_evaluate_Te)
        {
            this->field_evaluate_Te->evaluate(this->particle_group,
                                              Sym<REAL>("FLUID_TEMPERATURE"), 0,
                                              Te->GetCoeffs());
        }
        for (int d = 0; d < this->ndim; ++d)
        {
            if (this->field_evaluate_ve[d])
            {
                this->field_evaluate_ve[d]->evaluate(
                    this->particle_group, Sym<REAL>("FLUID_FLOW_SPEED"), d,
                    ve[d]->GetCoeffs());
            }
        }
    }

    inline void remove_marked_particles()
    {
        this->particle_remover->remove(this->particle_group,
                                       (*this->particle_group)[Sym<INT>("ID")],
                                       this->particle_remove_key);
    }

protected:
    virtual inline void integrate_inner_ion(ParticleSubGroupSharedPtr sg,
                                            const double dt_inner)
    {
        const auto k_dt = dt_inner;

        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:ions_3D", sg,
                [=](auto E, auto B, auto Q, auto M, auto P, auto V, auto TSP)
                {
                    const REAL dt_left  = k_dt - TSP.at(0);
                    const REAL hdt_left = dt_left * 0.5;
                    if (dt_left > 0.0)
                    {
                        const REAL QoM = Q.at(0) / M.at(0);

                        const REAL scaling_t = QoM * hdt_left;
                        const REAL t_0       = B.at(0) * scaling_t;
                        const REAL t_1       = B.at(1) * scaling_t;
                        const REAL t_2       = B.at(2) * scaling_t;

                        const REAL tmagsq = t_0 * t_0 + t_1 * t_1 + t_2 * t_2;
                        const REAL scaling_s = 2.0 / (1.0 + tmagsq);

                        const REAL s_0 = scaling_s * t_0;
                        const REAL s_1 = scaling_s * t_1;
                        const REAL s_2 = scaling_s * t_2;

                        const REAL V_0 = V.at(0);
                        const REAL V_1 = V.at(1);
                        const REAL V_2 = V.at(2);

                        const REAL v_minus_0 = V_0 + (E.at(0)) * scaling_t;
                        const REAL v_minus_1 = V_1 + (E.at(1)) * scaling_t;
                        const REAL v_minus_2 = V_2 + (E.at(2)) * scaling_t;

                        REAL v_prime_0, v_prime_1, v_prime_2;
                        MAPPING_CROSS_PRODUCT_3D(
                            v_minus_0, v_minus_1, v_minus_2, t_0, t_1, t_2,
                            v_prime_0, v_prime_1, v_prime_2)

                        v_prime_0 += v_minus_0;
                        v_prime_1 += v_minus_1;
                        v_prime_2 += v_minus_2;

                        REAL v_plus_0, v_plus_1, v_plus_2;
                        MAPPING_CROSS_PRODUCT_3D(v_prime_0, v_prime_1,
                                                 v_prime_2, s_0, s_1, s_2,
                                                 v_plus_0, v_plus_1, v_plus_2)

                        v_plus_0 += v_minus_0;
                        v_plus_1 += v_minus_1;
                        v_plus_2 += v_minus_2;

                        V.at(0) = v_plus_0 + scaling_t * (E.at(0));
                        V.at(1) = v_plus_1 + scaling_t * (E.at(1));
                        V.at(2) = v_plus_2 + scaling_t * (E.at(2));

                        // update of position to next time step
                        P.at(0) += dt_left * V.at(0);
                        P.at(1) += dt_left * V.at(1);
                        P.at(2) += dt_left * V.at(2);
                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("ELECTRIC_FIELD")),
                Access::read(Sym<REAL>("MAGNETIC_FIELD")),
                Access::read(Sym<REAL>("Q")), Access::read(Sym<REAL>("M")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
        else if (ndim == 2)
        {
            particle_loop(
                "ParticleSystem:ions_2D", sg,
                [=](auto E, auto B, auto Q, auto M, auto P, auto V, auto TSP)
                {
                    const REAL dt_left  = k_dt - TSP.at(0);
                    const REAL hdt_left = dt_left * 0.5;
                    if (dt_left > 0.0)
                    {
                        const REAL QoM = Q.at(0) / M.at(0);

                        const REAL scaling_t = QoM * hdt_left;
                        const REAL t_0       = B.at(0) * scaling_t;
                        const REAL t_1       = B.at(1) * scaling_t;
                        const REAL t_2       = B.at(2) * scaling_t;

                        const REAL tmagsq = t_0 * t_0 + t_1 * t_1 + t_2 * t_2;
                        const REAL scaling_s = 2.0 / (1.0 + tmagsq);

                        const REAL s_0 = scaling_s * t_0;
                        const REAL s_1 = scaling_s * t_1;
                        const REAL s_2 = scaling_s * t_2;

                        REAL o = hdt_left * V.at(2);
                        REAL h =
                            sycl::sqrt(1.0 + (o / P.at(0)) * (o / P.at(0)));

                        REAL V_0 = (V.at(0) + V.at(2) * o / P.at(0)) / h;
                        REAL V_1 = V.at(1);
                        REAL V_2 = (V.at(2) - V.at(0) * o / P.at(0)) / h;

                        const REAL v_minus_0 = V_0 + (E.at(0)) * scaling_t;
                        const REAL v_minus_1 = V_1 + (E.at(1)) * scaling_t;
                        const REAL v_minus_2 = V_2 + (E.at(2)) * scaling_t;

                        REAL v_prime_0, v_prime_1, v_prime_2;
                        MAPPING_CROSS_PRODUCT_3D(
                            v_minus_0, v_minus_1, v_minus_2, t_0, t_1, t_2,
                            v_prime_0, v_prime_1, v_prime_2)

                        v_prime_0 += v_minus_0;
                        v_prime_1 += v_minus_1;
                        v_prime_2 += v_minus_2;

                        REAL v_plus_0, v_plus_1, v_plus_2;
                        MAPPING_CROSS_PRODUCT_3D(v_prime_0, v_prime_1,
                                                 v_prime_2, s_0, s_1, s_2,
                                                 v_plus_0, v_plus_1, v_plus_2)

                        v_plus_0 += v_minus_0;
                        v_plus_1 += v_minus_1;
                        v_plus_2 += v_minus_2;

                        V_0 = v_plus_0 + scaling_t * (E.at(0));
                        V_1 = v_plus_1 + scaling_t * (E.at(1));
                        V_2 = v_plus_2 + scaling_t * (E.at(2));

                        // update of position to next time step
                        P.at(0) += dt_left * V_0;
                        P.at(1) += dt_left * V_1;

                        o = hdt_left * V_2;
                        h = sycl::sqrt(1.0 + (o / P.at(0)) * (o / P.at(0)));

                        V.at(0) = (V_0 + V_2 * o / P.at(0)) / h;
                        V.at(1) = V_1;
                        V.at(2) = (V_2 - V_0 * o / P.at(0)) / h;

                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("ELECTRIC_FIELD")),
                Access::read(Sym<REAL>("MAGNETIC_FIELD")),
                Access::read(Sym<REAL>("Q")), Access::read(Sym<REAL>("M")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
    }
    virtual inline void integrate_inner_neutral(ParticleSubGroupSharedPtr sg,
                                                const double dt_inner)
    {
        const auto k_dt = dt_inner;

        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:neutrals_3D", sg,
                [=](auto V, auto P, auto TSP)
                {
                    const REAL dt_left = k_dt - TSP.at(0);
                    if (dt_left > 0.0)
                    {
                        P.at(0) += dt_left * V.at(0);
                        P.at(1) += dt_left * V.at(1);
                        P.at(2) += dt_left * V.at(2);

                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
        else if (ndim == 2)
        {
            particle_loop(
                "ParticleSystem:neutrals_2D", sg,
                [=](auto V, auto P, auto TSP)
                {
                    const REAL dt_left  = k_dt - TSP.at(0);
                    const REAL hdt_left = dt_left * 0.5;

                    if (dt_left > 0.0)
                    {
                        REAL o = hdt_left * V.at(2);
                        REAL h =
                            sycl::sqrt(1.0 + (o / P.at(0)) * (o / P.at(0)));

                        REAL vx = (V.at(0) + V.at(2) * o / P.at(0)) / h;
                        REAL vz = (V.at(2) - V.at(0) * o / P.at(0)) / h;

                        P.at(0) += dt_left * vx;
                        P.at(1) += dt_left * V.at(1);

                        o = hdt_left * vz;
                        h = sycl::sqrt(1.0 + (o / P.at(0)) * (o / P.at(0)));

                        V.at(0) = (vx + vz * o / P.at(0)) / h;
                        V.at(2) = (vz - vx * o / P.at(0)) / h;

                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::write(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
    }

    virtual inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                const double dt_inner)
    {
        auto ions = particle_sub_group(
            sg, [=](auto Q) { return Q.at(0) != 0.0; },
            Access::read(Sym<REAL>("Q")));
        integrate_inner_ion(ions, dt_inner);
        auto neutrals = particle_sub_group(
            sg, [=](auto Q) { return Q.at(0) == 0.0; },
            Access::read(Sym<REAL>("Q")));
        integrate_inner_neutral(neutrals, dt_inner);
    }

    const long size;
    const long rank;
    std::mt19937 rng_phasespace;

    uint64_t total_num_particles_added = 0;

    const int particle_remove_key = -1;
    std::shared_ptr<ParticleRemover> particle_remover;

    std::shared_ptr<ParticleGroupTemporary> particle_group_temporary;

    std::map<std::string, SpeciesInfo> species_map;

    std::vector<Sym<REAL>> src_syms;
    std::vector<int> src_components;

    std::shared_ptr<FieldProject<DisContField>> field_project;

    std::shared_ptr<FieldProject<DisContField>> diagnostic_project;

    std::shared_ptr<FunctionEvaluateBasis<DisContField>> field_evaluate_ne;
    std::shared_ptr<FunctionEvaluateBasis<DisContField>> field_evaluate_Te;
    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_ve;

    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_E;

    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_B;

    std::shared_ptr<NektarCompositeTruncatedReflection> reflection;

    /// Simulation time
    double simulation_time;

    double mesh_length; // mesh conversion to m
    double Nnorm;       // Density normalisation to m^-3
    double Tnorm;       // Temperature normalisation to eV
    double Bnorm;       // B field normalisation to T
    double omega_c;     // Reference ion gyrofrequency

    inline void apply_timestep_reset()
    {
        particle_loop(
            this->particle_group,
            [=](auto TSP)
            {
                TSP.at(0) = 0.0;
                TSP.at(1) = 0.0;
            },
            Access::write(Sym<REAL>("TSP")))
            ->execute();
    }

    virtual void pre_advection(ParticleSubGroupSharedPtr sg)
    {
        reflection->pre_advection(sg);
    };

    virtual void apply_boundary_conditions(ParticleSubGroupSharedPtr sg,
                                           ParticleGroupSharedPtr cg, double dt)
    {
        reflection->execute(sg);
    };

    auto find_partial_moves(ParticleSubGroupSharedPtr sg, const double dt)
    {
        return particle_sub_group(
            sg, [=](auto TSP) { return TSP.at(0) < dt; },
            Access::read(Sym<REAL>("TSP")));
    };

    inline void apply_timestep_inner(const double dt)
    {
        auto child_group =
            this->particle_group_temporary->get(this->particle_group);
        auto sg = particle_sub_group(this->particle_group);
        pre_advection(sg);
        integrate_inner(sg, dt);
        apply_boundary_conditions(sg, child_group, dt);
        this->particle_group->add_particles_local(child_group);
        this->particle_group_temporary->restore(this->particle_group,
                                                child_group);
        sg = particle_sub_group(this->particle_group);
        sg = find_partial_moves(sg, dt);
        while (get_npart_global(sg) > 0)
        {
            auto child_group =
                this->particle_group_temporary->get(this->particle_group);
            pre_advection(sg);
            integrate_inner(sg, dt);
            apply_boundary_conditions(sg, child_group, dt);
            this->particle_group->add_particles_local(child_group);
            this->particle_group_temporary->restore(this->particle_group,
                                                    child_group);
            sg = particle_sub_group(this->particle_group);
            sg = find_partial_moves(sg, dt);
        }
    }

    inline void apply_timestep(const double dt)
    {
        apply_timestep_reset();

        apply_timestep_inner(dt);
    }

    /**
     *  Apply boundary conditions and transfer particles between MPI
     * ranks.
     * // Move some of this to PartSysBase / make it a pure-virtual
     * func?
     */
    inline void transfer_particles()
    {
        auto t0 = profile_timestamp();

        this->particle_group->hybrid_move();
        this->cell_id_translation->execute();
        this->particle_group->cell_move();
        this->sycl_target->profile_map.inc(
            "ParticleSystem", "transfer_particles", 1,
            profile_elapsed(t0, profile_timestamp()));
    }
};
} // namespace PENKNIFE
#endif
