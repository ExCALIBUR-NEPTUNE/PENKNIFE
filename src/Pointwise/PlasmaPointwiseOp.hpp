#pragma once

#include <MultiRegions/Common/Operator.hpp>

namespace PENKNIFE
{
using namespace Nektar;

// Pointwise equation-system calculations in the session's execution space.
// Fields retain their padded/interleaved blocks. Species stay inside the DoF
// loop.
template <typename TData>
class PlasmaPointwiseOp : public MultiRegions::Operator<TData>
{
public:
    using Field = LibUtilities::Field<TData, FieldState::Phys>;
    using Base  = MultiRegions::Operator<TData>;
    static inline const std::string name = "PlasmaPointwise";

    static std::shared_ptr<PlasmaPointwiseOp> Create(
        const MultiRegions::ExpListSharedPtr &expansion,
        const std::vector<std::string> &components,
        const std::string &exec = "")
    {
        return Base::template Create<PlasmaPointwiseOp>(expansion, components,
                                                        exec);
    }

    virtual void SetSources(const std::vector<unsigned> &density,
                            const std::vector<unsigned> &momentum,
                            const std::vector<unsigned> &energy,
                            unsigned electronEnergy)                 = 0;
    virtual void Quasineutrality(Field &in, Field &parallelCurrent, Field &ne,
                                 Field &Me)                          = 0;
    virtual void ParticleFields(Field &in, Field &ne, Field &Me, Field &b_unit,
                                Field &B, Field &E,
                                Field &m_particleFields)             = 0;
    virtual void AddSources(Field &out, Field &m_particleSources, Field &b_unit,
                            TData factor, TData timestep)            = 0;
    virtual void ParallelVelocities(Field &in, Field &b_unit, Field &ne,
                                    Field &Me, Field &m_velocity)    = 0;
    virtual void PressureWork(Field &in, Field &out, Field &m_velocityGradient,
                              Field &E, Field &b_unit, TData factor) = 0;
    virtual void PotentialForcing(Field &in, Field &b_unit, Field &mag_B,
                                  Field &m_phiRhs, Field &m_perpCoeff,
                                  unsigned omega_idx)                = 0;
    virtual void CorrectPotential(Field &in, Field &phi)             = 0;
    virtual void ElectricField(Field &m_phiGradient, Field &E)       = 0;
    virtual void DriftVelocities(Field &in, Field &m_gradient, Field &B,
                                 Field &mag_B, Field &E, Field &b_unit,
                                 Field &m_extraVelocity, Field &m_velocity,
                                 unsigned omega_idx)                 = 0;
    virtual void PolarizationTensor(Field &in, Field &m_velocity, Field &B,
                                    Field &mag_B, Field &m_tensor,
                                    TData omega_c)                   = 0;
    virtual void PolarizationFlux(Field &in, Field &m_tensorGradient,
                                  Field &m_omegaFlux)                = 0;
    virtual void CompleteVelocities(Field &in, Field &E, Field &B, Field &mag_B,
                                    Field &m_velocity)               = 0;

    void SetMagneticField(LibUtilities::Field<TData, FieldState::Phys> &b_unit)
    {
        m_b_unit = &b_unit;
    }
    void SetSpecies(const std::vector<unsigned> &density,
                    const std::vector<unsigned> &momentum,
                    const std::vector<unsigned> &energy,
                    const std::vector<TData> &masses,
                    const std::vector<TData> &charges)
    {
        v_SetSpecies(density, momentum, energy, masses, charges);
    }
    void SetElectronEnergyId(unsigned id)
    {
        v_SetElectronEnergyId(id);
    }

protected:
    virtual void v_SetSpecies(const std::vector<unsigned> &density,
                              const std::vector<unsigned> &momentum,
                              const std::vector<unsigned> &energy,
                              const std::vector<TData> &masses,
                              const std::vector<TData> &charges) = 0;
    virtual void v_SetElectronEnergyId(unsigned id)              = 0;
    PlasmaPointwiseOp(const Nektar::MultiRegions::ExpListSharedPtr &expansion,
                      const std::vector<std::string> &components)
        : Base(expansion, components)
    {
    }
    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit = nullptr;
};
} // namespace PENKNIFE
