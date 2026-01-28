#ifndef TOKAMAKSYSTEM_HPP
#define TOKAMAKSYSTEM_HPP

#include "../../NESO/include/nektar_interface/solver_base/neso_session_function.hpp"
#include "nektar_interface/solver_base/time_evolved_eqnsys_base.hpp"
#include "nektar_interface/utilities.hpp"

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>

#include "../BoundaryConditions/TokamakBndConds.hpp"
#include "../Misc/Constants.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"
#include "ImplicitHelper.hpp"
#include "MagneticField.hpp"

#include <solvers/solver_callback_handler.hpp>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;

namespace NESO::Solvers::tokamak
{

/**
 * @brief Equation system for the Tokamak proxyapp
 *
 */
class TokamakSystem
    : public TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>
{
    friend class MagneticField;
    friend class TokamakBaseBndCond;
    friend class VariableConverter;

public:
    friend class MemoryManager<TokamakSystem>;

    /// Name of class
    static std::string class_name;

    /**
     * @brief Create an instance of this class and initialise it.
     */
    static SU::EquationSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        SU::EquationSystemSharedPtr p =
            MemoryManager<TokamakSystem>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }

    virtual std::shared_ptr<ParticleSystem> GetParticleSystem();

    /// Callback handler to call user-defined callbacks
    SolverCallbackHandler<TokamakSystem> solver_callback_handler;

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
    TokamakSystem(const LU::SessionReaderSharedPtr &session,
                  const SD::MeshGraphSharedPtr &graph);

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

    Array<OneD, MR::ExpListSharedPtr> m_allfields;
    Array<OneD, MR::ExpListSharedPtr> m_saved;
    Array<OneD, MR::ExpListSharedPtr> m_indfields;
    int n_indep_fields;
    int n_species;
    int n_fields_per_species;

    /** Density source fields cast to DisContFieldSharedPtr for use in
     * particle evaluation/projection methods
     */
    std::vector<MR::DisContFieldSharedPtr> src_fields;
    std::vector<int> components;
    std::vector<Sym<REAL>> src_syms;

    /// Bool to enable/disable growth rate recordings
    bool energy_enstrophy_recording_enabled;

    VariableConverterSharedPtr m_varConv;

    /// Boundary Conditions
    std::shared_ptr<TokamakBoundaryConditions> m_bndConds;

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

    /// Number of particle timesteps per fluid timestep.
    int num_part_substeps;
    /// Number of time steps between particle trajectory step writes.
    int particle_output_freq;
    /// Particle timestep size.
    double part_timestep;

    std::shared_ptr<ImplicitHelper> m_implHelper;

    virtual void load_params() override;

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

} // namespace NESO::Solvers::tokamak
#endif
