#include "FluxVolAdvOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_FLUX_VOL_ADV(EXEC, TYPE)                                      \
    template <>                                                                \
    std::string FluxVolAdvOpImpl<NektarSpaces::EXEC, TYPE>::className =        \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            FluxVolAdvOp<TYPE>::name + #EXEC,                                  \
            FluxVolAdvOpImpl<NektarSpaces::EXEC, TYPE>::Instantiate,           \
            "Parallel advection")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_FLUX_VOL_ADV(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_ADV(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_ADV(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_FLUX_VOL_ADV(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_ADV(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_ADV(Device, float);
#endif
#endif

#undef REGISTER_FLUX_VOL_ADV_DIFF
} // namespace PENKNIFE::detail
