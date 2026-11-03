#include "FluxVolDiffOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_FLUX_VOL_DIFF(EXEC, TYPE)                                     \
    template <>                                                                \
    std::string FluxVolDiffOpImpl<                                             \
        NektarSpaces::EXEC, closures::Braginskii<TYPE>, TYPE>::className =     \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            FluxVolDiffOp<TYPE>::name + closures::Braginskii<TYPE>::name +     \
                #EXEC,                                                         \
            FluxVolDiffOpImpl<NektarSpaces::EXEC, closures::Braginskii<TYPE>,  \
                              TYPE>::Instantiate,                              \
            "Parallel Braginskii conduction and anomalous diffusion")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_FLUX_VOL_DIFF(Serial, double);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_DIFF(AVX, double);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_DIFF(Device, double);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_FLUX_VOL_DIFF(Serial, float);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_DIFF(AVX, float);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_DIFF(Device, float);
#endif
#endif

#undef REGISTER_FLUX_VOL_DIFF
} // namespace PENKNIFE::detail
