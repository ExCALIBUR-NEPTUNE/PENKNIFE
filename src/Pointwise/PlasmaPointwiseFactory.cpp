#include "PlasmaPointwiseOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_PLASMA_POINTWISE(EXEC, TYPE)                                  \
    template <>                                                                \
    std::string PlasmaPointwiseOpImpl<NektarSpaces::EXEC, TYPE>::className =   \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            PlasmaPointwiseOp<TYPE>::name + #EXEC,                             \
            PlasmaPointwiseOpImpl<NektarSpaces::EXEC, TYPE>::Instantiate,      \
            "Pointwise plasma equation-system calculations")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_PLASMA_POINTWISE(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_PLASMA_POINTWISE(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_PLASMA_POINTWISE(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_PLASMA_POINTWISE(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_PLASMA_POINTWISE(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_PLASMA_POINTWISE(Device, float);
#endif
#endif

#undef REGISTER_PLASMA_POINTWISE
} // namespace PENKNIFE::detail
