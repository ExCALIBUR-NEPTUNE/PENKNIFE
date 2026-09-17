#ifndef NULLSYSTEM_HPP
#define NULLSYSTEM_HPP
#include "PlasmaSystem.hpp"

namespace PENKNIFE
{
class NullSystem : public PlasmaSystem
{
public:
    friend class MemoryManager<NullSystem>;

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
            MemoryManager<NullSystem>::AllocateSharedPtr(session, graph);
        p->InitObject();
        return p;
    }
    ~NullSystem() override = default;

protected:
    NullSystem(const LU::SessionReaderSharedPtr &session,
               const SD::MeshGraphSharedPtr &graph);
    void load_params() override;
    void v_InitObject(bool DeclareFields = true) override;
    bool v_PostIntegrate(int step) override;

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);

    void DoParticles(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

private:
    int ee_idx;

    std::vector<int> ni_src_idx;
    std::vector<int> vi_src_idx;
    std::vector<int> ei_src_idx;

    std::map<std::string, Array<OneD, NekDouble>> ni;
    std::map<std::string, Array<OneD, NekDouble>> Ti;
    std::map<std::string, std::vector<Array<OneD, NekDouble>>> vi;

    NekDouble m_epsilon;
    bool m_useSpecVanVisc;
    NekDouble
        m_sVVCutoffRatio; // Cut-off ratio from which to start decaying modes
    NekDouble m_sVVDiffCoeff; // Diffusion coefficient of SVV modes
};

} // namespace PENKNIFE
#endif