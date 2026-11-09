#include "Closures/DiffusiveFields.hpp"
#include "FluxVolDiffOpImpl.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_FLUX_VOL_DIFF(EXEC, TYPE, MODEL)                              \
    template <>                                                                \
    std::string FluxVolDiffOpImpl<NektarSpaces::EXEC, closures::MODEL<TYPE>,   \
                                  TYPE>::className =                           \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            FluxVolDiffOp<TYPE>::name + closures::MODEL<TYPE>::name + #EXEC,   \
            FluxVolDiffOpImpl<NektarSpaces::EXEC, closures::MODEL<TYPE>,       \
                              TYPE>::Instantiate,                              \
            "Parallel Braginskii conduction and anomalous diffusion")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_FLUX_VOL_DIFF(Serial, double, Braginskii);
REGISTER_FLUX_VOL_DIFF(Serial, double, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(Serial, double, DoubleDiffusive);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_DIFF(AVX, double, Braginskii);
REGISTER_FLUX_VOL_DIFF(AVX, double, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(AVX, double, DoubleDiffusive);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_DIFF(Device, double, Braginskii);
REGISTER_FLUX_VOL_DIFF(Device, double, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(Device, double, DoubleDiffusive);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_FLUX_VOL_DIFF(Serial, float, Braginskii);
REGISTER_FLUX_VOL_DIFF(Serial, float, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(Serial, float, DoubleDiffusive);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_VOL_DIFF(AVX, float, Braginskii);
REGISTER_FLUX_VOL_DIFF(AVX, float, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(AVX, float, DoubleDiffusive);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_VOL_DIFF(Device, float, Braginskii);
REGISTER_FLUX_VOL_DIFF(Device, float, SingleDiffusive);
REGISTER_FLUX_VOL_DIFF(Device, float, DoubleDiffusive);
#endif
#endif

#undef REGISTER_FLUX_VOL_DIFF
} // namespace PENKNIFE::detail
