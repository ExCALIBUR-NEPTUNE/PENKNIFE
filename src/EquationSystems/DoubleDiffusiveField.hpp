#ifndef DOUBLEDIFFUSIVEFIELD_HPP
#define DOUBLEDIFFUSIVEFIELD_HPP
#include "../Misc/VariableConverter.hpp"
#include "PlasmaSystem.hpp"

namespace PENKNIFE
{
class DoubleDiffusiveField : public PlasmaSystem
{
public:
    friend class MemoryManager<DoubleDiffusiveField>;

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
            MemoryManager<DoubleDiffusiveField>::AllocateSharedPtr(session,
                                                                   graph);
        p->InitObject();
        return p;
    }
    ~DoubleDiffusiveField() override = default;

protected:
    DoubleDiffusiveField(const LU::SessionReaderSharedPtr &session,
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

    void CalcK(const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f);
    void CalcKappa(const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f);
    void CalcKappa(const Array<OneD, Array<OneD, NekDouble>> &in_arr);
    void CalcDiffTensor();

    void DoDiffusion(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

private:
    int pe_idx;

    std::vector<int> ni_src_idx;
    std::vector<int> ei_src_idx;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;
    NekDouble m_epsilon;

    bool m_useSpecVanVisc;
    NekDouble
        m_sVVCutoffRatio; // Cut-off ratio from which to start decaying modes
    NekDouble m_sVVDiffCoeff; // Diffusion coefficient of SVV modes

    NekDouble k_par;
    NekDouble k_perp;
    NekDouble k_cross;
    NekDouble kappa_i_par;
    NekDouble kappa_i_perp;
    NekDouble kappa_i_cross;
    NekDouble kappa_e_par;
    NekDouble kappa_e_perp;
    NekDouble kappa_e_cross;
    NekDouble k_ci;
    NekDouble k_ce;

    Array<OneD, NekDouble> m_kpar;
    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kcross;

    StdRegions::VarCoeffMap m_D;

    // For Diffusion
    // workaround for bug in DiffusionLDG
    Array<OneD, MR::ExpListSharedPtr> m_difffields;
};

} // namespace PENKNIFE
#endif