#ifndef REACTIONS_HPP
#define REACTIONS_HPP
#include "AMJUEL.hpp"
#include <neso_rng_toolkit.hpp>

using namespace VANTAGE::Reactions;

namespace PENKNIFE
{

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

namespace temp
{
// Temporary until component extraction is added to Reactions
struct ComponentDataOnDevice : public ReactionDataBaseOnDevice<1>
{

    ComponentDataOnDevice() = default;

    /**
     * @brief Function to extract particle dat values into an array
     *
     * @param index Read-only accessor to a loop index for a ParticleLoop
     * inside which calc_data is called. Access using either
     * index.get_loop_linear_index(), index.get_local_linear_index(),
     * index.get_sub_linear_index() as required.
     * @param req_int_props Vector of symbols for integer-valued properties that
     * need to be used for the reaction rate calculation.
     * @param req_real_props Vector of symbols for real-valued properties that
     * need to be used for the reaction rate calculation.
     * @param kernel The random number generator kernel potentially used in the
     * calculation
     *
     * @return A REAL-valued array of containing the extracted data
     */
    std::array<REAL, 1> calc_data(
        const Access::LoopIndex::Read &index,
        const Access::SymVector::Write<INT> &req_int_props,
        const Access::SymVector::Read<REAL> &req_real_props,
        typename ReactionDataBaseOnDevice<1>::RNG_KERNEL_TYPE::KernelType
            &kernel) const
    {

        std::array<REAL, 1> result;

        result[0] = req_real_props.at(this->prop_ind, index, this->comp_ind);

        return result;
    }

public:
    int prop_ind;
    int comp_ind;
};

/**
 * @brief Reaction data used to extract real valued ParticleDat
 */
struct ComponentData : public ReactionDataBase<ComponentDataOnDevice, 1>
{

    /**
     * @brief Constructor for ComponentData.
     *
     * @param extracted_sym The Sym<REAL> corresponding to the ParticleDat whose
     * components should be extracted
     * @param comp The component of the ParticleDat to be extracted
     */
    ComponentData(const Sym<REAL> &extracted_sym, const int comp)
        : ReactionDataBase<ComponentDataOnDevice, 1>(),
          extracted_sym(extracted_sym)
    {

        this->required_real_props.add(extracted_sym.name);
        this->on_device_obj = ComponentDataOnDevice();

        this->on_device_obj->comp_ind = comp;
        this->index_on_device_object();
    }

    /**
     * @brief Index the particle weight on the on-device object
     */
    void index_on_device_object()
    {

        this->on_device_obj->prop_ind =
            this->required_real_props.find_index(this->extracted_sym.name);
    };

private:
    Sym<REAL> extracted_sym;
};

auto inline component(const std::string &name, const int comp)
{

    return ComponentData(Sym<REAL>(name), comp);
}
} // namespace temp

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