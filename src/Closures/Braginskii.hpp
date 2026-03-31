#ifndef BRAGINSKII_HPP
#define BRAGINSKII_HPP

#include "Closure.hpp"

namespace SD = Nektar::SpatialDomains;

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
    void v_EvaluateClosure(
        const Array<OneD, Array<OneD, NekDouble>> &values,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &grads,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes,
        Array<OneD, Array<OneD, NekDouble>> &frictions,
        const Array<OneD, NekDouble> &ne,
        const Array<OneD, NekDouble> &ve) override;

    void CalcCollisionFrequencies(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, NekDouble> &ne);

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