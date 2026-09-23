#ifndef REACTIONS_HPP
#define REACTIONS_HPP
#include "AMJUEL.hpp"
#include <neso_rng_toolkit.hpp>

namespace PENKNIFE
{
using namespace VANTAGE::Reactions;
using namespace NP;

inline auto get_uniform_rng_kernel(SYCLTargetSharedPtr sycl_target,
                                   std::size_t n_samples,
                                   std::uint64_t root_seed = 141351)
{

    const int rank = sycl_target->comm_pair.rank_parent;

    std::uint64_t seed = NESO::RNGToolkit::create_seeds(
        sycl_target->comm_pair.size_parent, rank, root_seed);

    auto rng_uniform = NESO::RNGToolkit::create_rng<REAL>(
        NESO::RNGToolkit::Distribution::Uniform<REAL>{
            NESO::RNGToolkit::Distribution::next_value(0.0), 1.0},
        seed, sycl_target->device, sycl_target->device_index);

    // Create an interface between NESO-RNG-Toolkit and NESO-Particles KernelRNG
    auto rng_interface =
        make_rng_generation_function<GenericDeviceRNGGenerationFunction, REAL>(
            [=](REAL *d_ptr, const std::size_t num_samples) -> int
            { return rng_uniform->get_samples(d_ptr, num_samples); });

    auto rng_kernel =
        host_atomic_block_kernel_rng<REAL>(rng_interface, n_samples);

    return rng_kernel;
}

inline auto get_normal_rng_kernel(SYCLTargetSharedPtr sycl_target, REAL mean,
                                  REAL std, std::size_t n_samples,
                                  std::uint64_t root_seed = 1234561)
{

    const int rank = sycl_target->comm_pair.rank_parent;

    std::uint64_t seed = NESO::RNGToolkit::create_seeds(
        sycl_target->comm_pair.size_parent, rank, root_seed);

    auto rng_normal = NESO::RNGToolkit::create_rng<REAL>(
        NESO::RNGToolkit::Distribution::Normal<REAL>{mean, std}, seed,
        sycl_target->device, sycl_target->device_index);

    // Create an interface between NESO-RNG-Toolkit and NESO-Particles KernelRNG
    auto rng_interface =
        make_rng_generation_function<GenericDeviceRNGGenerationFunction, REAL>(
            [=](REAL *d_ptr, const std::size_t num_samples) -> int
            { return rng_normal->get_samples(d_ptr, num_samples); });

    auto rng_kernel =
        host_per_particle_block_rng<REAL>(rng_interface, n_samples);

    return rng_kernel;
}

template <size_t ndim, size_t vdim>
std::shared_ptr<AbstractReaction> specular_reflection(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate);

template <size_t vdim>
std::shared_ptr<AbstractReaction> surface_absorption(
    SYCLTargetSharedPtr sycl_target, Species absorbed_species, REAL rate);

template <size_t ndim, size_t vdim>
std::shared_ptr<AbstractReaction> thermal_reflection(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate,
    REAL std_dev);

template <size_t vdim>
std::shared_ptr<AbstractReaction> ionise_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target, double dens, double temp, double time,
    double vel, const Species &target_species, const Species &electron_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> ionise_reaction_fixed(
    SYCLTargetSharedPtr sycl_target, const Species &target_species,
    const Species &electron_species, REAL rate, REAL energy_rate);

template <size_t vdim>
std::shared_ptr<AbstractReaction> cx_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &projectile_species,
    const Species &target_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> cx_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &projectile_species, const Species &target_species, REAL rate,
    REAL sigma);

template <size_t vdim>
std::shared_ptr<AbstractReaction> recombination_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &marker_species,
    const Species &electron_species, const Species &neutral_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> recombination_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species, REAL rate, REAL energy_rate);

} // namespace PENKNIFE

#endif