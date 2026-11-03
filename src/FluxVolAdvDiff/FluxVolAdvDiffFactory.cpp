#include "FluxVolAdvDiffOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_FLUX_VOL_ADV_DIFF(EXEC, TYPE)                                 \
    template <>                                                                \
    std::string FluxVolAdvDiffOpImpl<                                          \
        NektarSpaces::EXEC, closures::Braginskii<TYPE>, TYPE>::className =     \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            FluxVolAdvDiffOp<TYPE>::name + closures::Braginskii<TYPE>::name +  \
                #EXEC,                                                         \
            FluxVolAdvDiffOpImpl<NektarSpaces::EXEC,                           \
                                 closures::Braginskii<TYPE>,                   \
                                 TYPE>::Instantiate,                           \
            "Combined parallel advection and Braginskii diffusion")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_FLUX_VOL_ADV_DIFF(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_ADV_DIFF(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_ADV_DIFF(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_FLUX_VOL_ADV_DIFF(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_ADV_DIFF(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_ADV_DIFF(Device, float);
#endif
#endif

#undef REGISTER_FLUX_VOL_ADV_DIFF
} // namespace PENKNIFE::detail
