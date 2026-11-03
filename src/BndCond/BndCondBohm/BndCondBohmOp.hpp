#pragma once

#include <boost/algorithm/string/case_conv.hpp>

#include <SolverCore/BndCond/BndCondUpdateOp.hpp>

namespace Nektar
{

template <typename TData>
class BndCondBohmOp : public SolverCore::BndCondUpdateOp<TData>
{
public:
    static std::shared_ptr<BndCondBohmOp<TData>> Create(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components,
        const std::string &execStr = "")
    {
        std::string EoSName;

        if (expansionList->GetSession()->DefinesEquationOfState())
        {
            EoSName = boost::to_upper_copy(
                expansionList->GetSession()->GetEquationOfState().type);
        }
        else
        {
            NEKERROR(ErrorUtil::efatal,
                     "No EquationOfState section defined in session file");
        }

        // The key carries the equation of state, as the volume flux operators
        // do, so the lookup cannot go through Operator::Create<TOperator>(),
        // which builds the key from the class name alone.
        const std::string execStr0 =
            (execStr == "") ? MultiRegions::Operator<TData>::GetOpExecSpace(
                                  expansionList->GetSession())
                            : execStr;

        const std::string requestedKey = name + EoSName + execStr0;

        MultiRegions::OperatorFactory<TData> &factory =
            MultiRegions::GetOperatorFactory<TData>();

        if (!factory.ModuleExists(requestedKey))
        {
            std::stringstream msg;
            msg << "No such operator: " << requestedKey << std::endl;
            factory.PrintAvailableClasses(msg);
            NEKERROR(ErrorUtil::efatal, msg.str());
        }

        return std::static_pointer_cast<BndCondBohmOp<TData>>(
            factory.CreateInstance(requestedKey, expansionList, components));
    }

    static inline const std::string name = "BndCondBohm";

protected:
    BndCondBohmOp(const MultiRegions::ExpListSharedPtr &expansionList,
                  const std::vector<std::string> &components)
        : SolverCore::BndCondUpdateOp<TData>(expansionList, components)
    {
    }

    ~BndCondBohmOp() override = default;
};

} // namespace Nektar
