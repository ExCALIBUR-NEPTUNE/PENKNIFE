#include "Reactions.hpp"

namespace PENKNIFE
{

template <size_t ndim, size_t vdim>
std::shared_ptr<AbstractReaction> specular_reflection(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate)
{

    auto rate_data = FixedRateData(rate);

    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto velocity_data       = ExtractorData<ndim>(Sym<REAL>("VELOCITY"));
    auto specular_reflection = SpecularReflectionData<ndim>();
    auto pipeline            = PipelineData(velocity_data, specular_reflection);
    auto reflection_kernels =
        LinearScatteringKernels<vdim>(reflected_species, properties_map);

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

template std::shared_ptr<AbstractReaction> specular_reflection<3, 3>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate);
template std::shared_ptr<AbstractReaction> specular_reflection<2, 3>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate);
template std::shared_ptr<AbstractReaction> specular_reflection<2, 2>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate);

template <size_t vdim>
std::shared_ptr<AbstractReaction> surface_absorption(
    SYCLTargetSharedPtr sycl_target, Species absorbed_species, REAL rate)
{
    auto rate_data = FixedRateData(rate);

    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.source_density] =
        "SURFACE_DENSITY_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto absorption_kernels =
        GeneralAbsorptionKernels<vdim>(absorbed_species, properties_map);

    auto absorption = std::make_shared<
        LinearReactionBase<0, FixedRateData, decltype(absorption_kernels)>>(
        sycl_target, absorbed_species.get_id(), std::array<int, 0>{}, rate_data,
        absorption_kernels);
    return absorption;
}

template std::shared_ptr<AbstractReaction> surface_absorption<3>(
    SYCLTargetSharedPtr sycl_target, Species absorbed_species, REAL rate);
template std::shared_ptr<AbstractReaction> surface_absorption<2>(
    SYCLTargetSharedPtr sycl_target, Species absorbed_species, REAL rate);

template <size_t ndim, size_t vdim>
std::shared_ptr<AbstractReaction> thermal_reflection(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate,
    REAL std_dev)
{
    auto rate_data = FixedRateData(rate);

    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.source_density] =
        "SURFACE_DENSITY_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_momentum] =
        "SURFACE_MOMENTUM_SOURCE";
    properties_map[VANTAGE::Reactions::default_properties.source_energy] =
        "SURFACE_ENERGY_SOURCE";

    auto cartesian_reflection = CartesianBasisReflectionData();

    auto reflection_kernels =
        LinearScatteringKernels<vdim>(reflected_species, properties_map);

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

template std::shared_ptr<AbstractReaction> thermal_reflection<3, 3>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate,
    REAL std_dev);
template std::shared_ptr<AbstractReaction> thermal_reflection<2, 3>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate,
    REAL std_dev);
template std::shared_ptr<AbstractReaction> thermal_reflection<2, 2>(
    SYCLTargetSharedPtr sycl_target, Species reflected_species, REAL rate,
    REAL std_dev);

template <size_t vdim>
std::shared_ptr<AbstractReaction> ionise_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target, double dens, double temp, double time,
    double vel, const Species &target_species, const Species &electron_species)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        "ELECTRON_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        "ELECTRON_TEMPERATURE";
    auto ionise_rate_data =
        AMJUEL::ionise_rate_data(dens, temp, time, properties_map);
    auto ionise_energy_data =
        AMJUEL::ionise_energy_data(dens, temp, time, vel, properties_map);

    auto ionise_reaction = std::make_shared<ElectronImpactIonisation<
        decltype(ionise_rate_data), decltype(ionise_energy_data), vdim>>(
        sycl_target, ionise_rate_data, ionise_energy_data, target_species,
        electron_species, properties_map);

    return ionise_reaction;
}

template std::shared_ptr<AbstractReaction> ionise_reaction_amjuel<3>(
    SYCLTargetSharedPtr sycl_target, double dens, double temp, double time,
    double vel, const Species &target_species, const Species &electron_species);

template std::shared_ptr<AbstractReaction> ionise_reaction_amjuel<2>(
    SYCLTargetSharedPtr sycl_target, double dens, double temp, double time,
    double vel, const Species &target_species, const Species &electron_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> ionise_reaction_fixed(
    SYCLTargetSharedPtr sycl_target, const Species &target_species,
    const Species &electron_species, REAL rate, REAL energy_rate)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        "ELECTRON_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        "ELECTRON_TEMPERATURE";
    auto ionise_rate_data   = FixedRateData(rate);
    auto ionise_energy_data = FixedRateData(energy_rate);

    auto ionise_reaction = std::make_shared<ElectronImpactIonisation<
        decltype(ionise_rate_data), decltype(ionise_energy_data), vdim>>(
        sycl_target, ionise_rate_data, ionise_energy_data, target_species,
        electron_species, properties_map);

    return ionise_reaction;
}

template std::shared_ptr<AbstractReaction> ionise_reaction_fixed<3>(
    SYCLTargetSharedPtr sycl_target, const Species &target_species,
    const Species &electron_species, REAL rate, REAL energy_rate);
template std::shared_ptr<AbstractReaction> ionise_reaction_fixed<2>(
    SYCLTargetSharedPtr sycl_target, const Species &target_species,
    const Species &electron_species, REAL rate, REAL energy_rate);

template <size_t vdim>
std::shared_ptr<AbstractReaction> cx_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &projectile_species,
    const Species &target_species)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        target_species.get_name() + "_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        target_species.get_name() + "_TEMPERATURE";
    properties_map[VANTAGE::Reactions::default_properties.fluid_flow_speed] =
        target_species.get_name() + "_FLOW_SPEED";

    auto parent_mass  = projectile_species.get_mass();
    auto child_mass   = target_species.get_mass();
    auto reduced_mass = (parent_mass * child_mass) / (parent_mass + child_mass);
    auto rate_data = AMJUEL::cx_rate_data(parent_mass, child_mass, dens, temp,
                                          time, vel, properties_map);
    auto cross_section = AMJUEL::amjuel_fit_cross_section(reduced_mass, vel);

    auto data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (child_mass * constants::mass_amu_SI * vel * vel),
            cross_section, rng_kernel, properties_map);

    auto data_calculator =
        DataCalculator<decltype(data_calc_sampler)>(data_calc_sampler);

    // The charge-exchange kernel, handles the descendant products and how
    // parent_species and descendant_species are modified by the reaction.
    auto cx_reaction_kernel = CXReactionKernels<vdim>(
        target_species, projectile_species, properties_map);

    // Designate that descendant particles have a "INTERNAL_STATE" that
    // corresponds to descendant_species
    const int out_state           = target_species.get_id();
    std::array<int, 1> out_states = {out_state};

    // Combining everything into a Reaction object
    auto cx_reaction = std::make_shared<
        LinearReactionBase<1, decltype(rate_data), decltype(cx_reaction_kernel),
                           decltype(data_calculator)>>(
        sycl_target, projectile_species.get_id(), out_states, rate_data,
        cx_reaction_kernel, data_calculator, properties_map);

    return cx_reaction;
}

template std::shared_ptr<AbstractReaction> cx_reaction_amjuel<3>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &projectile_species,
    const Species &target_species);
template std::shared_ptr<AbstractReaction> cx_reaction_amjuel<2>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &projectile_species,
    const Species &target_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> cx_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &projectile_species, const Species &target_species, REAL rate,
    REAL sigma)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        target_species.get_name() + "_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        target_species.get_name() + "_TEMPERATURE";
    properties_map[VANTAGE::Reactions::default_properties.fluid_flow_speed] =
        target_species.get_name() + "_FLOW_SPEED";

    auto rate_data     = FixedRateData(rate);
    auto cross_section = ConstantRateCrossSection(sigma);
    auto parent_mass   = projectile_species.get_mass();
    auto child_mass    = target_species.get_mass();
    auto reduced_mass = (parent_mass * child_mass) / (parent_mass + child_mass);

    auto data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (child_mass * constants::mass_amu_SI * vel * vel),
            cross_section, rng_kernel, properties_map);

    auto data_calculator =
        DataCalculator<decltype(data_calc_sampler)>(data_calc_sampler);

    // The charge-exchange kernel, handles the descendant products and how
    // parent_species and descendant_species are modified by the reaction.
    auto cx_reaction_kernel = CXReactionKernels<vdim>(
        target_species, projectile_species, properties_map);

    // Designate that descendant particles have a "INTERNAL_STATE" that
    // corresponds to descendant_species
    const int out_state           = target_species.get_id();
    std::array<int, 1> out_states = {out_state};

    // Combining everything into a Reaction object
    auto cx_reaction = std::make_shared<
        LinearReactionBase<1, decltype(rate_data), decltype(cx_reaction_kernel),
                           decltype(data_calculator)>>(
        sycl_target, projectile_species.get_id(), out_states, rate_data,
        cx_reaction_kernel, data_calculator, properties_map);

    return cx_reaction;
}

template std::shared_ptr<AbstractReaction> cx_reaction_fixed<3>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &projectile_species, const Species &target_species, REAL rate,
    REAL sigma);

template std::shared_ptr<AbstractReaction> cx_reaction_fixed<2>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &projectile_species, const Species &target_species, REAL rate,
    REAL sigma);

template <size_t vdim>
std::shared_ptr<AbstractReaction> recombination_reaction_amjuel(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &marker_species,
    const Species &electron_species, const Species &neutral_species)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        marker_species.get_name() + "_DENSITY";
    // properties_map[VANTAGE::Reactions::default_properties.fluid_temperature]
    // =
    //     marker_species.get_name() + "_TEMPERATURE";
    // properties_map[VANTAGE::Reactions::default_properties.fluid_flow_speed] =
    //     marker_species.get_name() + "_FLOW_SPEED";
    // properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
    //     "ELECTRON_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        "ELECTRON_TEMPERATURE";
    properties_map[VANTAGE::Reactions::default_properties.fluid_flow_speed] =
        "ELECTRON_FLOW_SPEED";

    auto recomb_data =
        AMJUEL::recomb_rate_data(dens, temp, time, properties_map);
    auto recomb_energy_data =
        AMJUEL::recomb_energy_data(dens, temp, time, vel, properties_map);

    auto constant_rate_cross_section = ConstantRateCrossSection(1.0);
    auto recomb_data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(constant_rate_cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (marker_species.get_mass() * constants::mass_amu_SI * vel *
                 vel),
            constant_rate_cross_section, rng_kernel, properties_map);
    auto recomb_data_calc_obj =
        DataCalculator<decltype(recomb_energy_data),
                       decltype(recomb_data_calc_sampler)>(
            recomb_energy_data, recomb_data_calc_sampler);

    double potential_energy =
        13.6 * constants::e / (constants::mass_amu_SI * vel * vel);

    auto recomb_reaction_kernel = RecombReactionKernels<vdim>(
        marker_species, electron_species, potential_energy, properties_map);

    const int out_state                  = neutral_species.get_id();
    std::array<int, 1> recomb_out_states = {out_state};

    auto recomb_reaction =
        std::make_shared<LinearReactionBase<1, decltype(recomb_data),
                                            decltype(recomb_reaction_kernel),
                                            decltype(recomb_data_calc_obj)>>(
            sycl_target, marker_species.get_id(), recomb_out_states,
            recomb_data, recomb_reaction_kernel, recomb_data_calc_obj,
            properties_map);
    return recomb_reaction;
}

template std::shared_ptr<AbstractReaction> recombination_reaction_amjuel<3>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &marker_species,
    const Species &electron_species, const Species &neutral_species);

template std::shared_ptr<AbstractReaction> recombination_reaction_amjuel<2>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double dens,
    double temp, double time, double vel, const Species &marker_species,
    const Species &electron_species, const Species &neutral_species);

template <size_t vdim>
std::shared_ptr<AbstractReaction> recombination_reaction_fixed(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species, REAL rate, REAL energy_rate)
{
    auto properties_map = get_default_map();
    properties_map[VANTAGE::Reactions::default_properties.fluid_density] =
        marker_species.get_name() + "_DENSITY";
    properties_map[VANTAGE::Reactions::default_properties.fluid_temperature] =
        marker_species.get_name() + "_TEMPERATURE";
    properties_map[VANTAGE::Reactions::default_properties.fluid_flow_speed] =
        marker_species.get_name() + "_FLOW_SPEED";

    auto recomb_data        = FixedRateData(rate);
    auto recomb_energy_data = FixedRateData(energy_rate);

    auto constant_rate_cross_section = ConstantRateCrossSection(1.0);
    auto recomb_data_calc_sampler =
        FilteredMaxwellianSampler<vdim, decltype(constant_rate_cross_section)>(
            (constants::temp_SI * constants::k_B) /
                (marker_species.get_mass() * constants::mass_amu_SI * vel *
                 vel),
            constant_rate_cross_section, rng_kernel, properties_map);
    auto recomb_data_calc_obj =
        DataCalculator<decltype(recomb_energy_data),
                       decltype(recomb_data_calc_sampler)>(
            recomb_energy_data, recomb_data_calc_sampler);

    double potential_energy =
        13.6 * constants::e / (constants::mass_amu_SI * vel * vel);

    auto recomb_reaction_kernel = RecombReactionKernels<vdim>(
        marker_species, electron_species, potential_energy, properties_map);

    const int out_state                  = neutral_species.get_id();
    std::array<int, 1> recomb_out_states = {out_state};

    auto recomb_reaction =
        std::make_shared<LinearReactionBase<1, decltype(recomb_data),
                                            decltype(recomb_reaction_kernel),
                                            decltype(recomb_data_calc_obj)>>(
            sycl_target, marker_species.get_id(), recomb_out_states,
            recomb_data, recomb_reaction_kernel, recomb_data_calc_obj,
            properties_map);
    return recomb_reaction;
}

template std::shared_ptr<AbstractReaction> recombination_reaction_fixed<3>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species, REAL rate, REAL energy_rate);

template std::shared_ptr<AbstractReaction> recombination_reaction_fixed<2>(
    SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<HostAtomicBlockKernelRNG<REAL>> rng_kernel, double vel,
    const Species &marker_species, const Species &electron_species,
    const Species &neutral_species, REAL rate, REAL energy_rate);

} // namespace PENKNIFE
