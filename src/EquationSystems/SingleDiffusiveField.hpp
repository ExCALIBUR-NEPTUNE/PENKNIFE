#ifndef SINGLEDIFFUSIVEFIELD_HPP
#define SINGLEDIFFUSIVEFIELD_HPP
#include "PlasmaSystem.hpp"

namespace PENKNIFE
{
class SingleDiffusiveField : public PlasmaSystem
{
public:
    friend class MemoryManager<SingleDiffusiveField>;

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
            MemoryManager<SingleDiffusiveField>::AllocateSharedPtr(session,
                                                                   graph);
        p->InitObject();
        return p;
    }
    ~SingleDiffusiveField() override = default;

protected:
    SingleDiffusiveField(const LU::SessionReaderSharedPtr &session,
                         const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool DeclareFields = true) override;
    bool v_PostIntegrate(int step) override;
    void ImplicitTimeIntCG(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);
    void CalcKPar(int f);
    void CalcKPerp(int f);
    void CalcKPerpAnomalous(int f);
    void CalcDiffTensor(int f);

    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

    // For Diffusion

    NekDouble m_epsilon;
    bool m_useSpecVanVisc;
    NekDouble
        m_sVVCutoffRatio; // Cut-off ratio from which to start decaying modes
    NekDouble m_sVVDiffCoeff; // Diffusion coefficient of SVV modes

    NekDouble k_par;
    NekDouble k_perp;
    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    Array<OneD, NekDouble> m_D[3][3];

    std::shared_ptr<DisContField> diag_field;
};

} // namespace PENKNIFE
#endif