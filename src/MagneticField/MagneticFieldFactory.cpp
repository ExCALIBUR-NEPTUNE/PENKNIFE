#include "MagneticFieldOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_MAGNETIC_FIELD(EXEC, TYPE)                                  \
    template <>                                                             \
    std::string MagneticFieldOpImpl<NektarSpaces::EXEC, TYPE>::className =     \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(    \
            MagneticFieldOp<TYPE>::name + #EXEC,                             \
            MagneticFieldOpImpl<NektarSpaces::EXEC, TYPE>::Instantiate,       \
            "Read a prescribed magnetic field and compute its unit vector")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_MAGNETIC_FIELD(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_MAGNETIC_FIELD(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_MAGNETIC_FIELD(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_MAGNETIC_FIELD(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_MAGNETIC_FIELD(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_MAGNETIC_FIELD(Device, float);
#endif
#endif

#undef REGISTER_MAGNETIC_FIELD
} // namespace PENKNIFE::detail
