#include "Closures/DiffusiveFields.hpp"
#include "FluxTraceAdvDiffOpImpl.hpp"
#include "RiemannSolvers/LaxFriedrichsSolverKernels.hpp"

namespace PENKNIFE::detail
{
// Match the execution spaces and precisions exported by the linked Nektar++.
#define REGISTER_FLUX_TRACE_ADV_DIFF(EXEC, TYPE, MODEL)                        \
    template <>                                                                \
    std::string FluxTraceAdvDiffOpImpl<LaxFriedrichsSolverKernel,              \
                                       closures::MODEL<TYPE>,                  \
                                       NektarSpaces::EXEC, TYPE>::className =  \
        MultiRegions::GetOperatorFactory<TYPE>().RegisterCreatorFunction(      \
            FluxTraceAdvDiffOp<TYPE>::name + std::string("LaxFriedrichs") +    \
                closures::MODEL<TYPE>::name + #EXEC,                           \
            FluxTraceAdvDiffOpImpl<LaxFriedrichsSolverKernel,                  \
                                   closures::MODEL<TYPE>, NektarSpaces::EXEC,  \
                                   TYPE>::Instantiate,                         \
            "Combined plasma advection and anisotropic trace diffusion")

#ifdef NEKTAR_ENABLE_DOUBLE_PRECISION
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, double, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, double, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, double, DoubleDiffusive);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, double, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, double, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, double, DoubleDiffusive);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_TRACE_ADV_DIFF(Device, double, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(Device, double, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(Device, double, DoubleDiffusive);
#endif
#endif

#ifdef NEKTAR_ENABLE_SINGLE_PRECISION
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, float, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, float, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(Serial, float, DoubleDiffusive);
#ifdef NEKTAR_ENABLE_SIMD
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, float, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, float, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(AVX, float, DoubleDiffusive);
#endif
#ifdef NEKTAR_ENABLE_DEVICE
REGISTER_FLUX_TRACE_ADV_DIFF(Device, float, Braginskii);
REGISTER_FLUX_TRACE_ADV_DIFF(Device, float, SingleDiffusive);
REGISTER_FLUX_TRACE_ADV_DIFF(Device, float, DoubleDiffusive);
#endif
#endif

#undef REGISTER_FLUX_TRACE_ADV_DIFF
} // namespace PENKNIFE::detail
