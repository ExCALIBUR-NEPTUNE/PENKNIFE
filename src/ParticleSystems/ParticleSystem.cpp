#include "ParticleSystem.hpp"

namespace NESO::Solvers::tokamak
{

std::string ParticleSystem::class_name =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ParticleSystem", ParticleSystem::create, "Particle System");

ParticleSystem::ParticleSystem(NESOReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph, MPI_Comm comm)
    : PartSysBase(session, graph, comm), simulation_time(0.0),
      size(this->sycl_target->comm_pair.size_parent),
      rank(this->sycl_target->comm_pair.rank_parent) {};

void ParticleSystem::set_up_species()
{
    // get seed from file
    std::srand(std::time(nullptr));
    int seed;

    this->config->get_session()->LoadParameter("particle_position_seed", seed,
                                               std::rand());
    this->rng_phasespace = std::mt19937(seed + this->rank);

    double particle_thermal_velocity;

    int s = 0;
    for (const auto &[k, v] : this->config->get_particle_species())
    {
        double particle_mass, particle_charge;
        this->config->load_particle_species_parameter(k, "Mass", particle_mass,
                                                      1.0);
        this->config->load_particle_species_parameter(k, "Charge",
                                                      particle_charge, 0.0);
        long particle_number = this->config->get_particle_species_initial_N(k);

        if (particle_number > 0)
        {
            std::vector<std::vector<double>> positions, velocities;
            std::vector<int> cells;
            auto vmap = this->config->get_particle_species_initial(k);
            if (auto v = vmap.find(std::pair("n", 0)); v != vmap.end())
            {
                rng_phasespace = dist_within_extents(
                    this->graph, v->second.m_expression, 0, particle_number,
                    positions, cells, 1.0e-10, this->rng_phasespace);
            }
            else // Point source
            {
                long rstart, rend;
                get_decomp_1d(this->size, particle_number, this->rank, &rstart,
                              &rend);
                const long local_particle_number = rend - rstart;

                double x = vmap.at(std::pair("X", 0)).m_expression->Evaluate();

                double y = vmap.at(std::pair("Y", 0)).m_expression->Evaluate();

                positions.emplace_back(
                    std::vector<double>(local_particle_number, x));

                positions.emplace_back(
                    std::vector<double>(local_particle_number, y));
                cells = std::vector<int>(local_particle_number, 0);

                if (ndim == 3)
                {
                    double z =
                        vmap.at(std::pair("Z", 0)).m_expression->Evaluate();
                    positions.emplace_back(
                        std::vector<double>(local_particle_number, z));
                }
            }

            int N         = cells.size();
            int id_offset = 0;
            MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM,
                              this->sycl_target->comm));
            if (N > 0)
            {
                double weight = 1.0;
                if (auto v = vmap.find(std::pair("W", 0)); v != vmap.end())
                {
                    weight = v->second.m_expression->Evaluate();
                }
                if (auto v = vmap.find(std::pair("T", 0)); v != vmap.end())
                {
                    double T   = v->second.m_expression->Evaluate();
                    double vth = constants::c *
                                 std::sqrt(Tnorm * T /
                                           (particle_mass * constants::m_p)) /
                                 (mesh_length * omega_c);
                    velocities = NESO::Particles::normal_distribution(
                        N, ndim, 0.0, vth, this->rng_phasespace);
                }

                else if (auto v = vmap.find(std::pair("Tin", 0));
                         v != vmap.end()) // Specific to EIRENE example
                {
                    for (int d = 0; d < 3; ++d)
                    {
                        velocities.emplace_back(std::vector<double>(N));
                    }

                    double T = v->second.m_expression->Evaluate();

                    std::uniform_real_distribution u(0.0, 1.0);
                    std::gamma_distribution mb(1.5, T);
                    std::normal_distribution norm(0.0,
                                                  std::sqrt(T / particle_mass));
                    double vth = constants::c *
                                 std::sqrt(Tnorm / constants::m_p) /
                                 (mesh_length * omega_c);

                    for (int p = 0; p < N; ++p)
                    {
                        double sintheta = std::sqrt(u(this->rng_phasespace));

                        // double energy = mb(this->rng_phasespace);
                        // double speed =
                        //     vth * std::sqrt(2 * energy / particle_mass);
                        // double phi       = 2 * M_PI *
                        // u(this->rng_phasespace); velocities[1][p] = speed *
                        // sintheta * cos(phi); velocities[2][p] = speed *
                        // sintheta * sin(phi);

                        // velocities[0][p] =
                        //     -speed * std::sqrt(1 - sintheta * sintheta);

                        velocities[1][p] = vth * norm(this->rng_phasespace);
                        velocities[2][p] = vth * norm(this->rng_phasespace);
                        double vperp =
                            std::sqrt(velocities[1][p] * velocities[1][p] +
                                      velocities[2][p] * velocities[2][p]);
                        velocities[0][p] = -vperp *
                                           std::sqrt(1 - sintheta * sintheta) /
                                           sintheta;
                    }
                }

                else // Explicit velocities
                {
                    if (auto v = vmap.find(std::pair("V", 0)); v != vmap.end())
                    {
                        particle_thermal_velocity =
                            v->second.m_expression->Evaluate();
                    }
                    if (auto v = vmap.find(std::pair("VX", 0)); v != vmap.end())
                    {
                        double vx = v->second.m_expression->Evaluate();
                        velocities.emplace_back(std::vector<double>(N, vx));
                    }
                    else
                    {
                        velocities.emplace_back(
                            NESO::Particles::normal_distribution(
                                N, 1, 0.0, particle_thermal_velocity,
                                this->rng_phasespace)[0]);
                    }
                    if (auto v = vmap.find(std::pair("VY", 0)); v != vmap.end())
                    {
                        double vy = v->second.m_expression->Evaluate();
                        velocities.emplace_back(std::vector<double>(N, vy));
                    }
                    else
                    {
                        velocities.emplace_back(
                            NESO::Particles::normal_distribution(
                                N, 1, 0.0, particle_thermal_velocity,
                                this->rng_phasespace)[0]);
                    }

                    if (auto v = vmap.find(std::pair("VZ", 0)); v != vmap.end())
                    {
                        double vz = v->second.m_expression->Evaluate();
                        velocities.emplace_back(std::vector<double>(N, vz));
                    }
                    else
                    {
                        velocities.emplace_back(
                            NESO::Particles::normal_distribution(
                                N, 1, 0.0, particle_thermal_velocity,
                                this->rng_phasespace)[0]);
                    }
                }
                ParticleSet initial_distribution(
                    N, this->particle_group->get_particle_spec());

                for (int px = 0; px < N; px++)
                {
                    for (int dimx = 0; dimx < this->graph->GetMeshDimension();
                         dimx++)
                    {
                        initial_distribution[Sym<REAL>("POSITION")][px][dimx] =
                            positions[dimx][px];
                        initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                            velocities[dimx][px];

                        initial_distribution[Sym<REAL>(
                            "ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
                        initial_distribution[Sym<REAL>("FLUID_FLOW_SPEED")][px]
                                            [dimx] = 0;
                    }
                    initial_distribution[Sym<REAL>("VELOCITY")][px][2] =
                        velocities[2][px];
                    initial_distribution[Sym<REAL>("Q")][px][0] =
                        particle_charge;
                    initial_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                    initial_distribution[Sym<INT>("ID")][px][0] =
                        px + id_offset + this->total_num_particles_added;
                    initial_distribution[Sym<INT>("CELL_ID")][px][0] =
                        cells.at(px);
                    initial_distribution[Sym<INT>("INTERNAL_STATE")][px][0] = s;
                    initial_distribution[Sym<REAL>("WEIGHT")][px][0] = weight;
                    initial_distribution[Sym<REAL>("TOT_REACTION_RATE")][px]
                                        [0] = 0.0;
                    initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                        2.0;
                    initial_distribution[Sym<REAL>("ELECTRON_TEMPERATURE")][px]
                                        [0] = 2.0;
                    initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")]
                                        [px][0] = 0.0;
                    initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")]
                                        [px][0] = 0.0;
                    initial_distribution[Sym<REAL>("FLUID_DENSITY")][px][0] =
                        2.0; // 1e18 m^-3
                    initial_distribution[Sym<REAL>("FLUID_TEMPERATURE")][px]
                                        [0] = 2.0; // eV
                }

                this->particle_group->add_particles_local(initial_distribution);
                this->total_num_particles_added += N;
            }
        }

        int state = s;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group,
                [state](auto sid) { return sid[0] == state; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        species_map[k] =
            SpeciesInfo{state, particle_mass, particle_charge, sub_group};
        s++;
    }
    set_up_boundaries();
}

template <typename RNG>
inline std::vector<double> gamma_distribution(const int N, const double alpha,
                                              const double beta, RNG &rng)
{
    std::gamma_distribution<> d{alpha, beta};
    std::vector<double> array(N);
    for (int px = 0; px < N; px++)
    {
        array[px] = d(rng);
    }

    return array;
}

void ParticleSystem::add_sources(double time, double dt)
{
    double particle_thermal_velocity;
    const long rank = this->sycl_target->comm_pair.rank_parent;

    for (auto &[k, v] : this->species_map)
    {
        double particle_mass   = v.mass;
        double particle_charge = v.charge;

        for (auto &source : this->config->get_particle_species_sources(k))
        {
            long particle_number = std::get<0>(source);
            auto &vmap           = std::get<3>(source);
            if (particle_number > 0)
            {
                double weight = 1.0;
                if (auto v = vmap.find(std::pair("W", 0)); v != vmap.end())
                {
                    weight = v->second.m_expression->Evaluate();
                }

                std::vector<std::vector<double>> positions, velocities;
                std::vector<int> cells;
                if (std::get<1>(source) ==
                    ParticleSourceType::eBulk) // Diffuse source
                {
                    auto v         = vmap.find(std::pair("n", 0));
                    rng_phasespace = dist_within_extents(
                        this->graph, v->second.m_expression, time,
                        particle_number, positions, cells, 1.0e-10,
                        this->rng_phasespace);
                }
                else if (std::get<1>(source) ==
                         ParticleSourceType::ePoint) // Point source
                {
                    long rstart, rend;
                    get_decomp_1d(this->size, particle_number, this->rank,
                                  &rstart, &rend);
                    const long local_particle_number = rend - rstart;

                    double x =
                        vmap.at(std::pair("X", 0)).m_expression->Evaluate();

                    double y =
                        vmap.at(std::pair("Y", 0)).m_expression->Evaluate();

                    positions.emplace_back(
                        std::vector<double>(local_particle_number, x));

                    positions.emplace_back(
                        std::vector<double>(local_particle_number, y));
                    cells = std::vector<int>(local_particle_number, 0);

                    if (ndim == 3)
                    {
                        double z =
                            vmap.at(std::pair("Z", 0)).m_expression->Evaluate();
                        positions.emplace_back(
                            std::vector<double>(local_particle_number, z));
                    }
                }
                else if (std::get<1>(source) == ParticleSourceType::eSurface)
                {
                    int br = *std::get<2>(source);
                }

                int N         = cells.size();
                int id_offset = 0;
                MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM,
                                  this->sycl_target->comm));
                if (N > 0)
                {
                    if (auto v = vmap.find(std::pair("T", 0)); v != vmap.end())
                    {
                        double T = v->second.m_expression->Evaluate();
                        double vth =
                            constants::c *
                            std::sqrt(Tnorm * T /
                                      (particle_mass * constants::m_p)) /
                            (mesh_length * omega_c);
                        velocities = NESO::Particles::normal_distribution(
                            N, ndim, 0.0, vth, this->rng_phasespace);
                    }

                    else if (auto v = vmap.find(std::pair("Tin", 0));
                             v != vmap.end()) // Specific to EIRENE example
                    {
                        for (int d = 0; d < 3; ++d)
                        {
                            velocities.emplace_back(std::vector<double>(N));
                        }
                        double T = v->second.m_expression->Evaluate();

                        std::uniform_real_distribution u(0.0, 1.0);
                        std::normal_distribution norm(0.0, 1.0);

                        double vth = constants::c *
                                     std::sqrt(Tnorm * T / (constants::m_p *
                                               particle_mass)) /
                                     (mesh_length * omega_c);

                        for (int p = 0; p < N; ++p)
                        {
                            velocities[1][p] = vth * norm(this->rng_phasespace);
                            velocities[2][p] = vth * norm(this->rng_phasespace);
                            velocities[0][p] =
                                -vth *
                                std::sqrt(-2 *
                                          std::log(u(this->rng_phasespace)));
                        }
                    }
                    else if (auto v = vmap.find(std::pair("Vin", 0));
                             v != vmap.end()) // Specific to EIRENE example
                    {
                        for (int d = 0; d < 3; ++d)
                        {
                            velocities.emplace_back(std::vector<double>(N));
                        }
                        double speed = v->second.m_expression->Evaluate()/
                                     (mesh_length * omega_c);;

                        std::uniform_real_distribution u(0.0, 1.0);

                        for (int p = 0; p < N; ++p)
                        {
                            // inverse transform sampling
                            double sintheta =
                                std::sqrt(u(this->rng_phasespace));
                            double phi = 2 * M_PI * u(this->rng_phasespace);
                            velocities[1][p] = speed * sintheta * cos(phi);
                            velocities[2][p] = speed * sintheta * sin(phi);

                            velocities[0][p] =
                                -speed * std::sqrt(1 - sintheta * sintheta);
                        }
                    }

                    else // Explicit velocities
                    {
                        if (auto v = vmap.find(std::pair("V", 0));
                            v != vmap.end())
                        {
                            particle_thermal_velocity =
                                v->second.m_expression->Evaluate();
                        }
                        if (auto v = vmap.find(std::pair("VX", 0));
                            v != vmap.end())
                        {
                            double vx = v->second.m_expression->Evaluate();
                            velocities.emplace_back(std::vector<double>(N, vx));
                        }
                        else
                        {
                            velocities.emplace_back(
                                NESO::Particles::normal_distribution(
                                    N, 1, 0.0, particle_thermal_velocity,
                                    this->rng_phasespace)[0]);
                        }
                        if (auto v = vmap.find(std::pair("VY", 0));
                            v != vmap.end())
                        {
                            double vy = v->second.m_expression->Evaluate();
                            velocities.emplace_back(std::vector<double>(N, vy));
                        }
                        else
                        {
                            velocities.emplace_back(
                                NESO::Particles::normal_distribution(
                                    N, 1, 0.0, particle_thermal_velocity,
                                    this->rng_phasespace)[0]);
                        }
                        if (this->ndim == 3)
                        {
                            if (auto v = vmap.find(std::pair("VZ", 0));
                                v != vmap.end())
                            {
                                double vz = v->second.m_expression->Evaluate();
                                velocities.emplace_back(
                                    std::vector<double>(N, vz));
                            }
                            else
                            {
                                velocities.emplace_back(
                                    NESO::Particles::normal_distribution(
                                        N, 1, 0.0, particle_thermal_velocity,
                                        this->rng_phasespace)[0]);
                            }
                        }
                    }
                    ParticleSet src_distribution(
                        N, this->particle_group->get_particle_spec());

                    for (int px = 0; px < N; px++)
                    {
                        for (int dimx = 0;
                             dimx < this->graph->GetMeshDimension(); dimx++)
                        {
                            src_distribution[Sym<REAL>("POSITION")][px][dimx] =
                                positions[dimx][px];
                            src_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                                velocities[dimx][px];

                            src_distribution[Sym<REAL>(
                                "ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
                            src_distribution[Sym<REAL>("FLUID_FLOW_SPEED")][px]
                                            [dimx] = 0;
                        }
                        src_distribution[Sym<REAL>("VELOCITY")][px][2] =
                            velocities[2][px];
                        src_distribution[Sym<REAL>("Q")][px][0] =
                            particle_charge;
                        src_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                        src_distribution[Sym<INT>("ID")][px][0] =
                            px + id_offset + this->total_num_particles_added;
                        src_distribution[Sym<INT>("CELL_ID")][px][0] =
                            cells.at(px);
                        src_distribution[Sym<INT>("INTERNAL_STATE")][px][0] =
                            v.id;
                        src_distribution[Sym<REAL>("WEIGHT")][px][0] = weight;
                        src_distribution[Sym<REAL>("TOT_REACTION_RATE")][px]
                                        [0] = 0.0;
                        src_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                            2.0;
                        src_distribution[Sym<REAL>("ELECTRON_TEMPERATURE")][px]
                                        [0] = 2.0;
                        src_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")]
                                        [px][0] = 0.0;
                        src_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")]
                                        [px][0] = 0.0;
                        src_distribution[Sym<REAL>("FLUID_DENSITY")][px][0] =
                            2.0; // 1e18 m^-3
                        src_distribution[Sym<REAL>("FLUID_TEMPERATURE")][px]
                                        [0] = 2.0; // eV
                    }

                    this->particle_group->add_particles_local(src_distribution);
                    this->total_num_particles_added += particle_number;
                }
            }
        }
        int s = v.id;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [s](auto sid) { return sid[0] == s; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        v.sub_group = sub_group;
    }
    transfer_particles();
}

void ParticleSystem::add_sinks(double time, double dt)
{
    std::uniform_real_distribution<> rng_dist(0, 1);

    Array<OneD, Array<OneD, NekDouble>> posarr(3);
    for (int dimx = 0; dimx < 3; dimx++)
    {
        posarr[dimx] = Array<OneD, NekDouble>(
            this->particle_group->get_npart_local(), 0.0);
    }

    // for each cell
    const int cell_count = this->domain->mesh->get_cell_count();
    for (int particle_index = 0, cellx = 0; cellx < cell_count; cellx++)
    {
        // for each particle in the cell

        auto cells =
            (*particle_group)[Sym<INT>("CELL_ID")]->cell_dat.get_cell(cellx);
        const int nrow = cells->nrow;
        // many cells will be empty so we check before issuing more
        // copies
        if (nrow > 0)
        {
            auto phys_positions =
                (*particle_group)[Sym<REAL>("POSITION")]->cell_dat.get_cell(
                    cellx);

            for (int rowx = 0; rowx < nrow; rowx++)
            {
                // copy the particle data into the store of points
                for (int dimx = 0; dimx < ndim; dimx++)
                {
                    posarr[dimx][particle_index] =
                        (*phys_positions)[dimx][rowx];
                }
                particle_index++;
            }
        }
    }
    int state = 0;
    for (const auto &[k, v] : this->config->get_particle_species())
    {
        for (auto &sink : this->config->get_particle_species_sinks(k))
        {
            LU::EquationSharedPtr neqn =
                sink.at(std::pair("n", 0)).m_expression;
            Array<OneD, NekDouble> eqnarr(
                this->particle_group->get_npart_local());

            neqn->Evaluate(posarr[0], posarr[1], posarr[2], time, eqnarr);

            for (int particle_index = 0, cellx = 0; cellx < cell_count; cellx++)
            {
                auto cells =
                    (*particle_group)[Sym<INT>("CELL_ID")]->cell_dat.get_cell(
                        cellx);
                const int nrow = cells->nrow;
                if (nrow > 0)
                {
                    auto id =
                        (*particle_group)[Sym<INT>("ID")]->cell_dat.get_cell(
                            cellx);
                    auto is = (*particle_group)[Sym<INT>("INTERNAL_STATE")]
                                  ->cell_dat.get_cell(cellx);
                    for (int rowx = 0; rowx < nrow; rowx++)
                    {
                        if (eqnarr[particle_index] >
                                rng_dist(this->rng_phasespace) &&
                            state == (*is)[0][rowx])
                        {
                            (*id)[0][rowx] = particle_remove_key;
                        }
                        particle_index++;
                    }
                    (*particle_group)[Sym<INT>("ID")]->cell_dat.set_cell(cellx,
                                                                         id);
                }
            }
        }
        state++;
    }
    remove_marked_particles();
}

void ParticleSystem::set_up_boundaries()
{
    auto store = std::make_shared<ParameterStore>();
    store->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                     1.0e-3);
    store->set<REAL>("CompositeIntersection/newton_tol", 1.0e-8);
    store->set<REAL>("CompositeIntersection/line_intersection_tol", 1.0e-10);
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);

    std::vector<int> reflection_composites;

    for (auto &[sk, sv] : this->config->get_particle_species_boundary(0))
    {
        if (sv == ParticleBoundaryConditionType::eReflective)
        {
            reflection_composites.push_back(sk);
        }
    }
    this->reflection = std::make_shared<NektarCompositeTruncatedReflection>(
        Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target, mesh,
        reflection_composites, store);
}

} // namespace NESO::Solvers::tokamak