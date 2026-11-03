#ifndef MAGNETIC_FIELD_OP
#define MAGNETIC_FIELD_OP

#include <MultiRegions/Common/Operator.hpp>

namespace PENKNIFE
{
using namespace Nektar;

template <typename TData>
class MagneticFieldOp : public MultiRegions::Operator<TData>
{
public:
    static std::shared_ptr<MagneticFieldOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        auto session = expansionList->GetSession();
        std::string execStr0 =
            (execStr == "")
                ? MultiRegions::Operator<TData>::GetOpExecSpace(session)
                : execStr;

        const std::string requestedKey = name + execStr0;

        MultiRegions::OperatorFactory<TData> &factory =
            MultiRegions::GetOperatorFactory<TData>();

        if (!factory.ModuleExists(requestedKey))
        {
            std::stringstream msg;
            msg << "No such operator: " << requestedKey << std::endl;
            factory.PrintAvailableClasses(msg);
            NEKERROR(ErrorUtil::efatal, msg.str());
        }

        return std::static_pointer_cast<MagneticFieldOp<TData>>(
            factory.CreateInstance(requestedKey, expansionList, components));
    }

    static inline const std::string name = "MagneticField";

    void Apply(const TData time = TData(0))
    {
        ASSERTL0(m_B && m_b_unit && m_mag_B,
                 "Bind B, unit B and squared magnitude with SetFields first.");
        v_Apply(time);
    }

    void operator()(const TData time = TData(0))
    {
        Apply(time);
    }

    // The equation system owns these fields and must outlive this operator.
    void SetFields(LibUtilities::Field<TData, FieldState::Phys> &B,
                   LibUtilities::Field<TData, FieldState::Phys> &b_unit,
                   LibUtilities::Field<TData, FieldState::Phys> &mag_B)
    {
        ASSERTL0(B.GetNumComponents() == 3 && b_unit.GetNumComponents() == 3 &&
                     mag_B.GetNumComponents() == 1,
                 "Magnetic fields require three vector components and one "
                 "squared-magnitude component.");
        ASSERTL0(B.GetNumHomoModes() == b_unit.GetNumHomoModes() &&
                     B.GetNumHomoModes() == mag_B.GetNumHomoModes(),
                 "Magnetic fields must have matching homogeneous modes.");
        ASSERTL0(B.GetBlocks().size() == b_unit.GetBlocks().size() &&
                     B.GetBlocks().size() == mag_B.GetBlocks().size(),
                 "Magnetic fields must share a block layout.");
        for (size_t blk = 0; blk < B.GetBlocks().size(); ++blk)
        {
            ASSERTL0(B.GetBlocks()[blk].CompSize() ==
                             b_unit.GetBlocks()[blk].CompSize() &&
                         B.GetBlocks()[blk].CompSize() ==
                             mag_B.GetBlocks()[blk].CompSize(),
                     "Magnetic fields must share a block layout.");
        }
        m_B      = &B;
        m_b_unit = &b_unit;
        m_mag_B  = &mag_B;
    }

protected:
    enum class field_type : int
    {
        mean_file,
        mean_exp,
        func,
        solve // Reserved for the later current-density solve.
    };

    field_type type                                        = field_type::func;
    LibUtilities::Field<TData, FieldState::Phys> *m_B      = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_b_unit = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_mag_B  = nullptr;

    MagneticFieldOp(const MultiRegions::ExpListSharedPtr &expansionList,
                    const std::vector<std::string> &components)
        : MultiRegions::Operator<TData>(expansionList, components)
    {
    }

    ~MagneticFieldOp() override = default;

    virtual void v_Apply(TData time) = 0;
};
} // namespace PENKNIFE
#endif
