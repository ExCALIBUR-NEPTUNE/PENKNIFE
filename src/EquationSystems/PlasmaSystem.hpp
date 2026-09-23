#ifndef PLASMA_SYSTEM_HPP
#define PLASMA_SYSTEM_HPP

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>

#include "../../NESO/include/nektar_interface/solver_base/neso_session_function.hpp"
#include "nektar_interface/utilities.hpp"
#include <solvers/solver_callback_handler.hpp>

#include "../BoundaryConditions/PlasmaBndConds.hpp"
#include "../Closures/Closure.hpp"
#include "../Misc/Constants.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"
#include "ImplicitHelper.hpp"
#include "MagneticField.hpp"



namespace PENKNIFE
{
using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;

/**
 * @brief Equation system for the PENKNIFE solver
 *
 */
class PlasmaSystem : public SU::UnsteadySystem
{
    friend class MagneticField;
    friend class PlasmaBaseBndCond;
    friend class VariableConverter;
    friend class Closure;

public:
    virtual std::shared_ptr<ParticleSystem> GetParticleSystem();

    /// Callback handler to call user-defined callbacks
    SolverCallbackHandler<PlasmaSystem> solver_callback_handler;

    struct Species
    {
        double mass;
        std::string name;
        std::map<int, int> fields;
    };
    std::map<int, Species> m_species;
    std::map<int, Species> &GetSpecies()
    {
        return m_species;
    }

    struct Ion : public Species
    {
        double charge;
    };
    std::map<int, Ion> m_ions;
    std::map<int, Ion> &GetIons()
    {
        return m_ions;
    }

    struct Neutral : public Species
    {
        int ion;
    };
    std::map<int, Neutral> m_neutrals;
    std::map<int, Neutral> &GetNeutrals()
    {
        return m_neutrals;
    }

protected:
    PlasmaSystem(const LU::SessionReaderSharedPtr &session,
                 const SD::MeshGraphSharedPtr &graph);

    NESOReaderSharedPtr neso_config;

    NESO::NektarFieldIndexMap field_to_index;

    /// Store mesh dims and number of quad points as member vars for convenience
    int n_dims;
    int n_pts;

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

    std::shared_ptr<MagneticField> mag_field;
    /// Magnetic field vector
    Array<OneD, MR::DisContFieldSharedPtr> B;
    /// Normalised magnetic field vector
    Array<OneD, Array<OneD, NekDouble>> b_unit;

    /// Squared Magnitude of the magnetic field
    Array<OneD, NekDouble> mag_B;

    /// Electric Field
    Array<OneD, MR::DisContFieldSharedPtr> E;

    /// Diffusion object used in anisotropic diffusion
    SU::DiffusionSharedPtr m_diffusion;

    MR::DisContFieldSharedPtr ne;
    MR::DisContFieldSharedPtr Te;
    Array<OneD, MR::DisContFieldSharedPtr> ve;

    Array<OneD, MR::ExpListSharedPtr> m_indfields;
    int n_indep_fields;
    int n_species;
    int n_fields_per_species;

    /** Density source fields cast to DisContFieldSharedPtr for use in
     * particle evaluation/projection methods
     */
    std::vector<MR::DisContFieldSharedPtr> src_fields;

    VariableConverterSharedPtr m_varConv;

    /// Boundary Conditions
    std::shared_ptr<PlasmaBoundaryConditions> m_bndConds;

    /// Forcing terms
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

    // Convenience key for varcoeffmaps
    static constexpr StdRegions::VarCoeffType vc[3][3] = {
        {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD01,
         StdRegions::eVarCoeffD02},
        {StdRegions::eVarCoeffD01, StdRegions::eVarCoeffD11,
         StdRegions::eVarCoeffD12},
        {StdRegions::eVarCoeffD02, StdRegions::eVarCoeffD12,
         StdRegions::eVarCoeffD22}};

    std::shared_ptr<ImplicitHelper> m_implHelper;

    virtual void load_params();

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    void SetBoundaryConditions(NekDouble time);

    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables) override;
    void v_DoSolve() override;
    virtual void v_DoInitialise(bool dump_initial_conditions) override;
    virtual void v_InitObject(bool DeclareField) override;
    virtual bool v_PostIntegrate(int step) override;
    virtual bool v_PreIntegrate(int step) override;
    virtual void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                        const int domain) override;

    NESOSessionFunctionSharedPtr get_species_function(
        const std::string &, std::string name,
        const MR::ExpListSharedPtr &field = MR::NullExpListSharedPtr,
        bool cache                        = false);
    std::map<std::string, std::map<std::string, NESOSessionFunctionSharedPtr>>
        m_nesoSessionFunctions;
};

} // namespace PENKNIFE
#endif
