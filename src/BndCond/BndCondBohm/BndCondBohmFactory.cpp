#include "BndCondBohmOpImpl.hpp"
namespace PENKNIFE::detail
{
#define REGISTER_BOHM(EXEC, TYPE)                                             \
    template <>                                                               \
    std::string BndCondBohmOpImpl<NektarSpaces::EXEC, TYPE>::className =      \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(     \
            BndCondBohmOp<TYPE>::name + #EXEC,                                \
            BndCondBohmOpImpl<NektarSpaces::EXEC, TYPE>::Instantiate,         \
            "Bohm-Chodura sheath boundary flux")
#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_BOHM(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_BOHM(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_BOHM(Device, double);
#endif
#endif
#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_BOHM(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_BOHM(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_BOHM(Device, float);
#endif
#endif
#undef REGISTER_BOHM
} // namespace PENKNIFE::detail
