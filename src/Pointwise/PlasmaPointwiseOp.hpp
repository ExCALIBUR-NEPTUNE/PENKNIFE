#pragma once
#include <MultiRegions/Common/Operator.hpp>
#include <limits>
namespace PENKNIFE
{
// Pointwise equation-system calculations in the session's execution space.
// Fields retain their padded/interleaved blocks. Species stay inside the DoF
// loop.
template <typename TData>
class PlasmaPointwiseOp : public Nektar::MultiRegions::Operator<TData>
{
public:
    using Field = Nektar::LibUtilities::Field<TData, Nektar::FieldState::Phys>;
    using Base  = Nektar::MultiRegions::Operator<TData>;
    static inline const std::string name = "PlasmaPointwise";
    static constexpr unsigned InvalidId  = std::numeric_limits<unsigned>::max();
    static std::shared_ptr<PlasmaPointwiseOp> Create(
        const Nektar::MultiRegions::ExpListSharedPtr &expansion,
        const std::vector<std::string> &components,
        const std::string &exec = "")
    {
        return Base::template Create<PlasmaPointwiseOp>(expansion, components,
                                                        exec);
    }
    virtual void SetSpecies(const std::vector<unsigned> &density,
                            const std::vector<unsigned> &momentum,
                            const std::vector<unsigned> &energy,
                            const std::vector<TData> &masses,
                            const std::vector<TData> &charges,
                            const std::vector<unsigned> &ions,
                            unsigned electronEnergy, TData electronMass,
                            TData backgroundTemperature)             = 0;
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

protected:
    PlasmaPointwiseOp(const Nektar::MultiRegions::ExpListSharedPtr &expansion,
                      const std::vector<std::string> &components)
        : Base(expansion, components)
    {
    }
};
} // namespace PENKNIFE
