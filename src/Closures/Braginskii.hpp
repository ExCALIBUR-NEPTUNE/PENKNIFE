#ifndef BRAGINSKII_HPP
#define BRAGINSKII_HPP

#include "Closure.hpp"

namespace PENKNIFE
{
// Forward declarations
class PlasmaSystem;
/**
 *
 */
class Braginskii : public Closure
{
public:
    friend class MemoryManager<Braginskii>;

    /// Creates an instance of this class
    static ClosureSharedPtr create(const std::weak_ptr<PlasmaSystem> &pSystem,
                                   const int spaceDim)
    {
        ClosureSharedPtr p =
            MemoryManager<Braginskii>::AllocateSharedPtr(pSystem, spaceDim);
        return p;
    }

    /// Name of the class
    static std::string className;

private:
    void v_CollisionFrequencies(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, NekDouble> &ne) override;

    void v_EvaluateHeatFlux(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes) override;

    void v_EvaluateThermalForce(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        const Array<OneD, NekDouble> &ne,
        Array<OneD, Array<OneD, NekDouble>> &force) override;

    void v_EvaluateFrictionHeating(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, NekDouble> &ne, const Array<OneD, NekDouble> &ve,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        Array<OneD, Array<OneD, NekDouble>> &heats) override;

    void v_EvaluateConductivity(const Array<OneD, NekDouble> &ne,
                                Array<OneD, NekDouble> &sigma) override;

    Braginskii(const std::weak_ptr<PlasmaSystem> &pSystem, const int spaceDim);

    ~Braginskii() = default;

    std::map<std::pair<int, int>, Array<OneD, NekDouble>> nu_ii;
    std::map<int, Array<OneD, NekDouble>> nu_ei;
    Array<OneD, NekDouble> nu_ee;

    std::map<int, Array<OneD, NekDouble>> nu_i;
    Array<OneD, NekDouble> nu_e;

    NekDouble k_ci;
    NekDouble k_ce;
};

} // namespace PENKNIFE

#endif