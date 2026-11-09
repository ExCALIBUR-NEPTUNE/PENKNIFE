#pragma once
#include <SolverCore/Advection/AdvectionVolumeFluxOp.hpp>
#include <limits>

namespace PENKNIFE
{
using namespace Nektar;
template <typename T>
class FluxVolAdvOp : public SolverCore::AdvectionVolumeFluxOp<T>
{
public:
    using Field = LibUtilities::Field<T, FieldState::Phys>;
    static std::shared_ptr<FluxVolAdvOp> Create(
        const MultiRegions::ExpListSharedPtr &exp,
        const std::vector<std::string> &names, const std::string & = "")
    {
        return std::make_shared<FluxVolAdvOp>(exp, names);
    }
    FluxVolAdvOp(const MultiRegions::ExpListSharedPtr &exp,
                 const std::vector<std::string> &names)
        : SolverCore::AdvectionVolumeFluxOp<T>(exp, names),
          dim(exp->GetCoordim(0))
    {
    }
    void SetMagneticField(Field &field)
    {
        magnetic = &field;
    }
    void SetElectronFields(Field &density, Field &momentum, T mass,
                           unsigned energy)
    {
        ne = &density;
        Me = &momentum;
        me = mass;
        ee = energy;
    }
    void SetIndices(const std::vector<unsigned> &n,
                    const std::vector<unsigned> &v,
                    const std::vector<unsigned> &e, const std::vector<T> &m,
                    const std::vector<T> &)
    {
        nIds   = n;
        vIds   = v;
        eIds   = e;
        masses = m;
    }
    void SetExtraFlux(Field &velocity, Field &flux, unsigned omega)
    {
        extraVelocity = &velocity;
        omegaFlux     = &flux;
        omegaId       = omega;
    }

protected:
    void v_Apply(Field &in, Field &out) override
    {
        out.SetInterleaveWidth(in);
        for (auto *f : {magnetic, ne, Me, extraVelocity, omegaFlux})
            if (f)
                f->SetInterleaveWidth(in);
        for (unsigned blk = 0; blk < in.GetBlocks().size(); ++blk)
        {
            const size_t stride = in.GetBlocks()[blk].CompSize();
            auto u              = in.GetBlocks()[blk]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto b = magnetic->GetBlocks()[blk]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto n = ne->GetBlocks()[blk]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto M = Me->GetBlocks()[blk]
                         .template GetPtr<NektarSpaces::HostSpace, ReadOnly>();
            auto dst =
                this->m_append
                    ? out.GetBlocks()[blk]
                          .template GetPtr<NektarSpaces::HostSpace, ReadWrite>()
                    : out.GetBlocks()[blk]
                          .template GetPtr<NektarSpaces::HostSpace,
                                           WriteOnly>();
            const T *extra =
                extraVelocity
                    ? extraVelocity->GetBlocks()[blk]
                          .template GetPtr<NektarSpaces::HostSpace, ReadOnly>()
                    : nullptr;
            const T *omega =
                omegaFlux
                    ? omegaFlux->GetBlocks()[blk]
                          .template GetPtr<NektarSpaces::HostSpace, ReadOnly>()
                    : nullptr;
            if (!this->m_append)
                std::fill(dst, dst + stride * out.GetNumComponents(), T(0));
            for (size_t p = 0; p < stride; ++p)
            {
                for (unsigned s = 0; s < nIds.size(); ++s)
                {
                    if (nIds[s] == invalid)
                        continue;
                    T density = u[nIds[s] * stride + p];
                    T momentum =
                        vIds[s] == invalid ? 0 : u[vIds[s] * stride + p];
                    T energy = eIds[s] == invalid ? 0 : u[eIds[s] * stride + p];
                    T velocity =
                        density > 0 ? momentum / (masses[s] * density) : 0;
                    for (unsigned d = 0; d < dim; ++d)
                    {
                        T ud = velocity * b[d * stride + p];
                        dst[(nIds[s] * dim + d) * stride + p] +=
                            this->m_scale * density * ud;
                        if (vIds[s] != invalid)
                            dst[(vIds[s] * dim + d) * stride + p] +=
                                this->m_scale *
                                (momentum * ud +
                                 (T(2) / 3) * energy * b[d * stride + p]);
                        if (eIds[s] != invalid)
                            dst[(eIds[s] * dim + d) * stride + p] +=
                                this->m_scale * energy * ud;
                    }
                }
                T ue = n[p] > 0 ? M[p] / (me * n[p]) : 0;
                for (unsigned d = 0; d < dim; ++d)
                    dst[(ee * dim + d) * stride + p] += this->m_scale *
                                                        u[ee * stride + p] *
                                                        ue * b[d * stride + p];
                if (extra)
                    for (unsigned f = 0; f < in.GetNumComponents(); ++f)
                        for (unsigned d = 0; d < dim; ++d)
                            dst[(f * dim + d) * stride + p] +=
                                this->m_scale *
                                (f == omegaId
                                     ? omega[d * stride + p]
                                     : u[f * stride + p] *
                                           extra[(f * dim + d) * stride + p]);
            }
        }
    }

private:
    static constexpr unsigned invalid = std::numeric_limits<unsigned>::max();
    unsigned dim, ee = invalid, omegaId = invalid;
    T me            = 1;
    Field *magnetic = nullptr, *ne = nullptr, *Me = nullptr,
          *extraVelocity = nullptr, *omegaFlux = nullptr;
    std::vector<unsigned> nIds, vIds, eIds;
    std::vector<T> masses;
};
} // namespace PENKNIFE
