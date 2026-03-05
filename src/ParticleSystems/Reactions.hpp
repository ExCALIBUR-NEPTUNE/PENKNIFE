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
inline auto specular_reflection(SYCLTargetSharedPtr sycl_target,
                                Species reflected_species, REAL rate)
{

    auto rate_data = FixedRateData(rate);

    auto properties_map = PropertiesMap();
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto velocity_data       = ExtractorData<ndim>(Sym<REAL>("VELOCITY"));
    auto specular_reflection = SpecularReflectionData<ndim>();
    auto pipeline            = PipelineData(velocity_data, specular_reflection);
    auto reflection_kernels  = LinearScatteringKernels<vdim>(
        reflected_species, properties_map.get_map());

    if constexpr (ndim == 2 && vdim == 3)
    {
        auto velocity_data_2 = temp::ComponentData(Sym<REAL>("VELOCITY"), 2);
        auto concat          = ConcatenatorData(pipeline, velocity_data_2);
        auto data_calculator = DataCalculator<decltype(concat)>(concat);

        auto reflection = std::make_shared<
            LinearReactionBase<1, FixedRateData, decltype(reflection_kernels),
                               decltype(data_calculator)>>(
            sycl_target, reflected_species.get_id(),
            std::array<int, 1>{static_cast<int>(reflected_species.get_id())},
            rate_data, reflection_kernels, data_calculator);
        return reflection;
    }
    else
    {
        auto data_calculator = DataCalculator<decltype(pipeline)>(pipeline);

        auto reflection = std::make_shared<
            LinearReactionBase<1, FixedRateData, decltype(reflection_kernels),
                               decltype(data_calculator)>>(
            sycl_target, reflected_species.get_id(),
            std::array<int, 1>{static_cast<int>(reflected_species.get_id())},
            rate_data, reflection_kernels, data_calculator);
        return reflection;
    }
}

template <size_t vdim>
inline auto surface_absorption(SYCLTargetSharedPtr sycl_target,
                               Species absorbed_species, REAL rate)
{
    auto rate_data = FixedRateData(rate);

    auto properties_map = PropertiesMap();
    properties_map[VANTAGE::Reactions::default_properties.source_density] =
        "SURFACE_DENSITY_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto absorption_kernels = GeneralAbsorptionKernels<vdim>(
        absorbed_species, properties_map.get_map());

    auto absorption = std::make_shared<
        LinearReactionBase<0, FixedRateData, decltype(absorption_kernels)>>(
        sycl_target, absorbed_species.get_id(), std::array<int, 0>{}, rate_data,
        absorption_kernels);
    return absorption;
}

template <size_t ndim, size_t vdim>
inline auto thermal_reflection(SYCLTargetSharedPtr sycl_target,
                               Species reflected_species, REAL rate,
                               REAL std_dev)
{
    auto rate_data = FixedRateData(rate);

    auto properties_map = PropertiesMap();
    properties_map[VANTAGE::Reactions::default_properties.source_density] =
        "SURFACE_DENSITY_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto cartesian_reflection = CartesianBasisReflectionData();

    auto reflection_kernels = LinearScatteringKernels<vdim>(
        reflected_species, properties_map.get_map());

    auto sampler1 =
        SamplerData(get_normal_rng_kernel(sycl_target, 0, std_dev, 1, 123456));
    auto sampler2 =
        SamplerData(get_normal_rng_kernel(sycl_target, 0, std_dev, 1, 654321));

    auto sampler_uniform =
        SamplerData(get_uniform_rng_kernel(sycl_target, 1, 987654));

    auto unary_lambda = [=](const REAL &U)
    { return std_dev * Kernel::sqrt(-2 * Kernel::log(U)); };

    auto lambda_wrapper = utils::LambdaWrapper(unary_lambda);
    auto unary_transform_data =
        uetData<1, decltype(lambda_wrapper)>(lambda_wrapper);
    auto rayleigh_sample = PipelineData(sampler_uniform, unary_transform_data);
    auto velocities     = ConcatenatorData(sampler1, sampler2, rayleigh_sample);
    auto reflected_data = PipelineData(velocities, cartesian_reflection);
    auto data_calculator = DataCalculator(reflected_data);
    auto reflection      = std::make_shared<
             LinearReactionBase<1, FixedRateData, decltype(reflection_kernels),
                                decltype(data_calculator)>>(
        sycl_target, reflected_species.get_id(),
        std::array<int, 1>{static_cast<int>(reflected_species.get_id())},
        rate_data, reflection_kernels, data_calculator);

    return reflection;
}

template <size_t vdim>
inline auto ionise_reaction_amjuel(SYCLTargetSharedPtr sycl_target,
                                   const Species &target_species,
                                   const Species &electron_species)
{
    auto ionise_rate_data   = AMJUEL::ionise_rate_data();
    auto ionise_energy_data = AMJUEL::ionise_energy_data();

    auto ionise_reaction = std::make_shared<ElectronImpactIonisation<
        decltype(ionise_rate_data), decltype(ionise_energy_data), vdim>>(
        sycl_target, ionise_rate_data, ionise_energy_data, target_species,
        electron_species);

    return ionise_reaction;
}

template <size_t vdim>
inline auto ionise_reaction_fixed(SYCLTargetSharedPtr sycl_target,
                                  const Species &target_species,
                                  const Species &electron_species, REAL rate,
                                  REAL energy_rate)
{
    auto ionise_rate_data   = FixedRateData(rate);
    auto ionise_energy_data = FixedRateData(energy_rate);
    auto ionise_reaction    = std::make_shared<ElectronImpactIonisation<
           decltype(ionise_rate_data), decltype(ionise_energy_data), vdim>>(
        sycl_target, ionise_rate_data, ionise_energy_data, target_species,
        electron_species);

    return ionise_reaction;
}

template <size_t vdim>
inline auto cx_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel,
    const Species &parent_species, const Species &descendant_species)
{
    auto parent_mass  = parent_species.get_mass();
    auto child_mass   = descendant_species.get_mass();
    auto reduced_mass = (parent_mass * child_mass) / (parent_mass + child_mass);
    auto rate_data    = AMJUEL::cx_rate_data(parent_mass, child_mass);
    auto cross_section = AMJUEL::amjuel_fit_cross_section(reduced_mass);

    auto data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (child_mass * norm::mass_amu_SI * norm::vel * norm::vel),
            cross_section, rng_kernel);

    auto data_calculator =
        DataCalculator<decltype(data_calc_sampler)>(data_calc_sampler);

    // The charge-exchange kernel, handles the descendant products and how
    // parent_species and descendant_species are modified by the reaction.
    auto cx_reaction_kernel =
        CXReactionKernels<vdim>(descendant_species, parent_species);

    // Designate that descendant particles have a "INTERNAL_STATE" that
    // corresponds to descendant_species
    const int out_state           = descendant_species.get_id();
    std::array<int, 1> out_states = {out_state};

    // Combining everything into a Reaction object
    auto cx_reaction = std::make_shared<
        LinearReactionBase<1, decltype(rate_data), decltype(cx_reaction_kernel),
                           decltype(data_calculator)>>(
        sycl_target, parent_species.get_id(), out_states, rate_data,
        cx_reaction_kernel, data_calculator);

    return cx_reaction;
}

template <size_t vdim>
inline auto cx_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel,
    const Species &parent_species, const Species &descendant_species, REAL rate,
    REAL sigma)
{
    auto rate_data     = FixedRateData(rate);
    auto cross_section = ConstantRateCrossSection(sigma);
    auto parent_mass   = parent_species.get_mass();
    auto child_mass    = descendant_species.get_mass();
    auto reduced_mass = (parent_mass * child_mass) / (parent_mass + child_mass);

    auto data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (child_mass * norm::mass_amu_SI * norm::vel * norm::vel),
            cross_section, rng_kernel);

    auto data_calculator =
        DataCalculator<decltype(data_calc_sampler)>(data_calc_sampler);

    // The charge-exchange kernel, handles the descendant products and how
    // parent_species and descendant_species are modified by the reaction.
    auto cx_reaction_kernel =
        CXReactionKernels<vdim>(descendant_species, parent_species);

    // Designate that descendant particles have a "INTERNAL_STATE" that
    // corresponds to descendant_species
    const int out_state           = descendant_species.get_id();
    std::array<int, 1> out_states = {out_state};

    // Combining everything into a Reaction object
    auto cx_reaction = std::make_shared<
        LinearReactionBase<1, decltype(rate_data), decltype(cx_reaction_kernel),
                           decltype(data_calculator)>>(
        sycl_target, parent_species.get_id(), out_states, rate_data,
        cx_reaction_kernel, data_calculator);

    return cx_reaction;
}

template <size_t vdim>
inline auto recombination_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species)
{
    auto recomb_data        = AMJUEL::recomb_rate_data();
    auto recomb_energy_data = AMJUEL::recomb_energy_data();

    auto constant_rate_cross_section = ConstantRateCrossSection(1.0);
    auto recomb_data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(constant_rate_cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (marker_species.get_mass() * norm::mass_amu_SI * norm::vel *
                 norm::vel),
            constant_rate_cross_section, rng_kernel);
    auto recomb_data_calc_obj =
        DataCalculator<decltype(recomb_energy_data),
                       decltype(recomb_data_calc_sampler)>(
            recomb_energy_data, recomb_data_calc_sampler);

    auto recomb_reaction_kernel = RecombReactionKernels<vdim>(
        marker_species, electron_species, norm::potential_energy);

    const int out_state                  = neutral_species.get_id();
    std::array<int, 1> recomb_out_states = {out_state};

    auto recomb_reaction =
        std::make_shared<LinearReactionBase<1, decltype(recomb_data),
                                            decltype(recomb_reaction_kernel),
                                            decltype(recomb_data_calc_obj)>>(
            sycl_target, marker_species.get_id(), recomb_out_states,
            recomb_data, recomb_reaction_kernel, recomb_data_calc_obj);
    return recomb_reaction;
}

template <size_t vdim>
inline auto recombination_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species, REAL rate, REAL energy_rate)
{
    auto recomb_data        = FixedRateData(rate);
    auto recomb_energy_data = FixedRateData(energy_rate);

    auto constant_rate_cross_section = ConstantRateCrossSection(1.0);
    auto recomb_data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(constant_rate_cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (marker_species.get_mass() * norm::mass_amu_SI * norm::vel *
                 norm::vel),
            constant_rate_cross_section, rng_kernel);
    auto recomb_data_calc_obj =
        DataCalculator<decltype(recomb_energy_data),
                       decltype(recomb_data_calc_sampler)>(
            recomb_energy_data, recomb_data_calc_sampler);

    auto recomb_reaction_kernel = RecombReactionKernels<vdim>(
        marker_species, electron_species, norm::potential_energy);

    const int out_state                  = neutral_species.get_id();
    std::array<int, 1> recomb_out_states = {out_state};

    auto recomb_reaction =
        std::make_shared<LinearReactionBase<1, decltype(recomb_data),
                                            decltype(recomb_reaction_kernel),
                                            decltype(recomb_data_calc_obj)>>(
            sycl_target, marker_species.get_id(), recomb_out_states,
            recomb_data, recomb_reaction_kernel, recomb_data_calc_obj);
    return recomb_reaction;
}

} // namespace PENKNIFE

#endif