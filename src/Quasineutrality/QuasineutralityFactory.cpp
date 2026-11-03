#include "QuasineutralityOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_QUASINEUTRALITY(EXEC, TYPE)                                   \
    template <>                                                                \
    std::string QuasineutralityOpImpl<NektarSpaces::EXEC, TYPE>::className =   \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            QuasineutralityOp<TYPE>::name + #EXEC,                             \
            QuasineutralityOpImpl<NektarSpaces::EXEC, TYPE>::Instantiate,      \
            "Calculate electron density by quasineutrality")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_QUASINEUTRALITY(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_QUASINEUTRALITY(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_QUASINEUTRALITY(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_QUASINEUTRALITY(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_QUASINEUTRALITY(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_QUASINEUTRALITY(Device, float);
#endif
#endif

#undef REGISTER_QUASINEUTRALITY
} // namespace PENKNIFE::detail
