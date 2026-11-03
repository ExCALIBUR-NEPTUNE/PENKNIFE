#ifndef PLASMA_SYSTEM_HPP
#define PLASMA_SYSTEM_HPP

#include "../../NESO/include/nektar_interface/solver_base/neso_session_function.hpp"
#include "nektar_interface/utilities.hpp"

#include "../BndCond/PlasmaBndConds.hpp"
#include "../MagneticField/MagneticFieldOp.hpp"
#include "../Quasineutrality/QuasineutralityOp.hpp"
#include "../Misc/Constants.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"
#include "ImplicitHelper.hpp"

#include <limits>
#include <vector>

#include <LibUtilities/BasicUtils/Field/Field.hpp>
#include <MultiRegions/ElmtOps/Expression/ExpressionOp.hpp>
#include <MultiRegions/ElmtOps/PhysTraceExtract/PhysTraceExtractOp.hpp>
#include <SolverCore/EquationSystems/UnsteadySystem.h>
#include <SolverCore/Forcing/Forcing.h>

namespace PENKNIFE
{
using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SC = Nektar::SolverCore;

/**
 * @brief Equation system for the PENKNIFE solver
 *
 */
class PlasmaSystem : public SC::UnsteadySystem
{
    friend class MagneticField;
    friend class PlasmaBaseBndCond;
    friend class VariableConverter;
    friend class Closure;

public:
    virtual std::shared_ptr<ParticleSystem> GetParticleSystem();

    // One global species index addresses every metadata vector. Ion and
    // neutral lists contain these indices, not positions within a subset.
    static constexpr unsigned int InvalidId =
        std::numeric_limits<unsigned int>::max();

    unsigned int GetNumSpecies() const
    {
        return m_names.size();
    }
    const std::vector<std::string> &GetSpeciesNames() const
    {
        return m_names;
    }
    const std::vector<unsigned int> &GetIons() const
    {
        return m_ions;
    }
    const std::vector<unsigned int> &GetNeutrals() const
    {
        return m_neutrals;
    }
    const std::vector<double> &GetMasses() const
    {
        return m_masses;
    }
    const std::vector<double> &GetCharges() const
    {
        return m_charges;
    }
    const std::vector<unsigned int> &GetDensityIds() const
    {
        return n_ids;
    }
    const std::vector<unsigned int> &GetMomentumIds() const
    {
        return v_ids;
    }
    const std::vector<unsigned int> &GetEnergyIds() const
    {
        return e_ids;
    }
    bool HasDensity(unsigned int s) const
    {
        return n_ids.at(s) != InvalidId;
    }
    bool HasMomentum(unsigned int s) const
    {
        return v_ids.at(s) != InvalidId;
    }
    bool HasEnergy(unsigned int s) const
    {
        return e_ids.at(s) != InvalidId;
    }
    unsigned int GetDensityId(unsigned int s) const
    {
        ASSERTL0(HasDensity(s),
                 "Species " + m_names.at(s) + " has no density field.");
        return n_ids[s];
    }
    unsigned int GetMomentumId(unsigned int s) const
    {
        ASSERTL0(HasMomentum(s),
                 "Species " + m_names.at(s) + " has no momentum field.");
        return v_ids[s];
    }
    unsigned int GetEnergyId(unsigned int s) const
    {
        ASSERTL0(HasEnergy(s),
                 "Species " + m_names.at(s) + " has no energy field.");
        return e_ids[s];
    }
    unsigned int GetFieldId(unsigned int s, const std::string &field) const
    {
        if (field == "n")
            return GetDensityId(s);
        if (field == "v")
            return GetMomentumId(s);
        if (field == "e")
            return GetEnergyId(s);
        ASSERTL0(false,
                 "Unknown species field '" + field + "'; expected n, v or e.");
        return InvalidId;
    }
    unsigned int GetNeutralIonId(unsigned int s) const
    {
        ASSERTL0(m_neutral_ions.at(s) != InvalidId,
                 "Species " + m_names.at(s) + " has no associated ion.");
        return m_neutral_ions[s];
    }

    void Quasineutrality(LU::Field<double, FieldState::Phys> &in);

protected:
    PlasmaSystem(const LU::SessionReaderSharedPtr &session,
                 const SD::MeshGraphSharedPtr &graph);

    NESOReaderSharedPtr neso_config;

    NESO::NektarFieldIndexMap field_to_index;

    /// Particle system
    std::shared_ptr<ParticleSystem> particle_sys;

    /// Flag identifying whether particles were enabled in the config file
    bool particles_enabled;

    /// List of field names required by the solver
    std::vector<std::string> required_fld_names;

    NekDouble mesh_length; // mesh conversion to m
    NekDouble Nnorm;       // Density normalisation to m^-3
    NekDouble Tnorm;       // Temperature normalisation to eV
    NekDouble Bnorm;       // B field normalisation to T
    NekDouble omega_c;     // Reference ion gyrofrequency
    NekDouble rho_s;       // Reference ion length scale
    NekDouble cs;          // Reference ion sound speed
    NekDouble me;          // Electron mass

    std::map<int, std::map<int, CompositeSharedPtr>> m_domains;
    std::map<int, int> m_dom_to_offset;

    bool transient_field;

    std::shared_ptr<MagneticFieldOp<double>> mag_field;

    /// Magnetic field vector
    LU::Field<double, FieldState::Phys> B;
    /// Normalised magnetic field vector
    LU::Field<double, FieldState::Phys> b_unit;
    LU::Field<double, FieldState::Phys> b_unit_trace;
    std::shared_ptr<PhysTraceExtractOp<double>> m_BTraceExtractOp;

    /// Squared Magnitude of the magnetic field
    LU::Field<double, FieldState::Phys> mag_B;

    /// Electric Field
    LU::Field<double, FieldState::Phys> E;

    LU::Field<double, FieldState::Phys> ne;
    /// Parallel current normalized by e Nnorm cs; prescribed zero for now.
    LU::Field<double, FieldState::Phys> j_par;
    /// Parallel electron momentum density: M_e = m_e n_e u_e_parallel.
    LU::Field<double, FieldState::Phys> Me;
    std::shared_ptr<QuasineutralityOp<double>> m_quasineutralityOp;

    unsigned int n_species = 0;
    std::vector<std::string> m_names;
    std::vector<double> m_masses;
    std::vector<double> m_charges;
    std::vector<unsigned int> m_ions;
    std::vector<unsigned int> m_neutrals;
    // Indexed by global species ID; InvalidId for species without an ion link.
    std::vector<unsigned int> m_neutral_ions;

    // Exactly one entry per species, including InvalidId for absent fields.
    std::vector<unsigned int> n_ids;
    std::vector<unsigned int> v_ids;
    std::vector<unsigned int> e_ids;
    int ee_idx;

    /** Density source fields cast to DisContFieldSharedPtr for use in
     * particle evaluation/projection methods
     */
    LU::Field<double, FieldState::Phys> src_fields;

    /// Forcing terms
    std::vector<SC::ForcingSharedPtr> m_forcing;

    // Convenience key for varcoeffmaps
    static constexpr StdRegions::VarCoeffType vc[3][3] = {
        {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD01,
         StdRegions::eVarCoeffD02},
        {StdRegions::eVarCoeffD01, StdRegions::eVarCoeffD11,
         StdRegions::eVarCoeffD12},
        {StdRegions::eVarCoeffD02, StdRegions::eVarCoeffD12,
         StdRegions::eVarCoeffD22}};

    std::shared_ptr<ImplicitHelper> m_implHelper;

    virtual void LoadParams();

    virtual void v_GenerateSummary(SummaryList &s) override;
    virtual void v_WriteFld(const std::string &) override;
    virtual void v_DoInitialise(bool dump_initial_conditions) override;

    virtual void v_InitObject(bool declareExpansionLists = true) override;
    virtual void v_InitialiseFields() override;
    virtual void v_InitialiseOperators() override;

    virtual bool v_PostIntegrate() override;
    virtual bool v_PreIntegrate() override;
    virtual void v_SetInitialConditions(NekDouble init_time) override;

    NESOSessionFunctionSharedPtr get_species_function(
        const std::string &, std::string name,
        const MR::ExpListSharedPtr &field = MR::NullExpListSharedPtr,
        bool cache                        = false);
    std::map<std::string, std::map<std::string, NESOSessionFunctionSharedPtr>>
        m_nesoSessionFunctions;
}; // namespace Nektar::SolverCore

} // namespace PENKNIFE
#endif
