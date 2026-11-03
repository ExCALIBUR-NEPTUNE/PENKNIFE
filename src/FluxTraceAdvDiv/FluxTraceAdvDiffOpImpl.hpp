#pragma once
#include "FluxTraceAdvDiffOp.hpp"

#include "FluxTraceAdvDiv/FluxTraceAdvKernels.hpp"
#include "FluxTraceAdvDiv/FluxTraceDiffKernels.hpp"
#include "LibUtilities/BasicUtils/Math/MathKernels.hpp"
#include "SolverCore/RiemannSolver/RiemannSolverKernels.hpp"
#include "SolverCore/TraceFlux/TraceFluxOpImpl.hpp"

namespace PENKNIFE::detail
{

/**
 * Combined plasma advection and anisotropic diffusion trace flux.
 * The inviscid launcher supplies F_adv; the diffusion kernel produces
 * -F_adv + F_diff.n, using the volume closure at the arithmetic mean state.
 * Interior/partition traces support either trace formulation. Boundaries
 * retain the diffusion ghost state and the boundary gradient path.
 */
template <template <typename> typename RiemannKernel, typename ClosureType,
          typename ExecSpace, typename TData>
class FluxTraceAdvDiffOpImpl
    : public SolverCore::detail::TraceFluxOpImpl<ExecSpace, TData,
                                                 FluxTraceAdvDiffOp<TData>>
{
    using MemSpace = typename ExecSpace::memory_space;
    using Base = SolverCore::detail::TraceFluxOpImpl<ExecSpace, TData,
                                                     FluxTraceAdvDiffOp<TData>>;

  public:
    /**
     * @brief Construct the operator, build the equation of state, size the
     * workspaces and read the transport parameters.
     *
     * @param expansionList - Expansion list the operator acts on.
     * @param components    - Names of the conserved variables.
     */
    FluxTraceAdvDiffOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                           const std::vector<std::string> &components)
        : SolverCore::detail::TraceFluxOpImpl<ExecSpace, TData,
                                              FluxTraceAdvDiffOp<TData>>(
              expansionList, components)
    {
        m_nDim = expansionList->GetCoordim(0);
        ASSERTL0(
            m_nDim == expansionList->GetExp(0)->GetShapeDimension(),
            "Embedded manifolds are not supported by this trace operator.");
        m_nComp = components.size();
        LibUtilities::SessionReaderSharedPtr session =
            expansionList->GetSession();

        session->LoadParameter("Nnorm", m_parameters.Nnorm, 1e18);
        session->LoadParameter("Tnorm", m_parameters.Tnorm, 100.);
        session->LoadParameter("mesh_length", m_parameters.mesh_length, 1.);
        TData Bnorm;
        session->LoadParameter("Bnorm", Bnorm, 1.);
        m_parameters.omega_c = TData(constants::qeomp) * Bnorm;
        for (auto value : {m_parameters.Nnorm, m_parameters.Tnorm,
                           m_parameters.mesh_length, m_parameters.omega_c})
            ASSERTL0(std::isfinite(value) && value > 0,
                     "Diffusion normalization must be finite and positive.");
        session->LoadParameter("Anomalous_D", m_anomD, 0.);
        session->LoadParameter("Anomalous_nu", m_anomNu, 0.);
        session->LoadParameter("Anomalous_chi", m_anomChi, 0.);
        session->LoadParameter("IPPenaltyCoeff", m_penaltyCoeff, 1.);
        ASSERTL0(std::isfinite(m_penaltyCoeff) && m_penaltyCoeff >= 0,
                 "IPPenaltyCoeff must be finite and nonnegative.");

        size_t maxNTraceNPts = 0;
        for (unsigned b = 0; b < this->m_intT0.size(); ++b)
        {
            maxNTraceNPts =
                std::max(maxNTraceNPts, this->m_intT0[b].m_nTraceXnPtsPad);
        }

        for (unsigned b = 0; b < this->m_bndT0.size(); ++b)
        {
            maxNTraceNPts =
                std::max(maxNTraceNPts, this->m_bndT0[b].m_nTraceXnPtsPad);
        }

        for (unsigned b = 0; b < this->m_parT0.size(); ++b)
        {
            maxNTraceNPts =
                std::max(maxNTraceNPts, this->m_parT0[b].m_nTraceXnPtsPad);
        }
        // Runtime-sized workspace: no fixed limit on the species count.
        const size_t stateSize =
            this->PadToVectorWidth(maxNTraceNPts * m_nComp);
        const size_t gradientSize =
            this->PadToVectorWidth(maxNTraceNPts * m_nComp * m_nDim);
        m_wspDiffusion =
            LibUtilities::MemoryRegion<TData>(stateSize + 2 * gradientSize);
        m_wspDiffusion.template Initialize<MemSpace>(0);
        auto *diffBase = m_wspDiffusion.template GetPtr<MemSpace, ReadWrite>(
            this->m_opStreamID);
        m_diffusionWorkspace = {diffBase, diffBase + stateSize,
                                diffBase + stateSize + gradientSize};

        // One length per array, summed to size the buffer and then walked to
        // carve the pointers, so the total cannot fall behind the carving. A
        // chunk appended to the carving alone puts the whole of it past the
        // end of the buffer, where it neither faults nor holds its values -
        // it reads back as whatever else owns that memory, which is a long
        // way from the NaNs it eventually produces. The assertion is the
        // cheap way to keep that from being a debugging exercise, so it is
        // ASSERTL0: it runs once per construction and the failure it catches
        // is silent in a release build.
        const size_t nDerivPts = maxNTraceNPts * m_nDim * m_nComp;
        const size_t nScalarPts = maxNTraceNPts;

        // MemoryRegion rather than std::vector: these arrays reach the AVX
        // kernels, which reinterpret_cast them to tinysimd::simd<TData> and
        // use aligned moves. std::vector only guarantees max_align_t (16 bytes
        // here), which is short of the 32 AVX2 wants and the 64 AVX512 wants,
        // so its base can fault before any offset is even added. MemoryRegion
        // defaults to NektarSpaces::host_memory_alignment for this reason.
        // The strides are padded so each array also starts on a vector.
        const size_t derivStride = this->PadToVectorWidth(nDerivPts);
        const size_t scalarStride = this->PadToVectorWidth(nScalarPts);

        m_wspDeriv = LibUtilities::MemoryRegion<TData>(
            2 * derivStride + scalarStride +
            this->PadToVectorWidth(maxNTraceNPts * m_nComp));
        m_wspDeriv.template Initialize<MemSpace>(0);

        TData *base = m_wspDeriv.template GetPtr<MemSpace, WriteOnly>();
        TData *cursor = base;
        m_gloDerivT0 = cursor;
        cursor += derivStride;
        m_gloDerivT1 = cursor;
        cursor += derivStride;

        m_penFactor = cursor;
        cursor += scalarStride;
        m_bndFluxWeight = cursor;
        cursor += this->PadToVectorWidth(maxNTraceNPts * m_nComp);

        ASSERTL0(cursor == base + m_wspDeriv.size(),
                 "m_wspDeriv is not the exact sum of the arrays carved out of "
                 "it - add the new array's length to the size above.");

        AllocateAuxiliary(this->m_intT0, m_intAux);
        AllocateAuxiliary(this->m_parT0, m_parAux);
        AllocateAuxiliary(this->m_bndT0, m_bndAux);
        // The exterior state the diffusion sees, which is not always the one
        // the Riemann solver sees - see the boundary loop in v_Apply().
        m_wspGloT1Diff = LibUtilities::MemoryRegion<TData>(
            this->PadToVectorWidth(maxNTraceNPts * m_nComp));
        m_wspGloT1Diff.template Initialize<MemSpace>(0);
        m_gloT1Diff = m_wspGloT1Diff.template GetPtr<MemSpace, WriteOnly>();
    }

    /**
     * @brief Ask each boundary condition what it does to the boundary flux.
     *
     * What kind of exterior state it built, and what weight it puts on the
     * diffusive flux. Both are the condition's own business - this operator
     * only supplies the trace layout to apply them in - so both are asked
     * rather than inferred from the boundary tag.
     *
     * Done on the host once: the mapping from trace to region is fixed after
     * construction, and so is the set of attached conditions.
     */
    void BuildBndFluxWeights()
    {
        m_bndFluxWeightPerTrace.assign(this->m_bndTraceIdT1.size(), {});
        m_bndGhostIsReflected.assign(this->m_bndTraceIdT1.size(), {});
        m_bndCondDiffState.assign(this->m_bndTraceIdT1.size(), {});
        m_bndCondDiffStateBlk.assign(this->m_bndTraceIdT1.size(), {});

        for (unsigned b = 0; b < this->m_bndTraceIdT1.size(); ++b)
        {
            const auto numTrace = this->m_bndTraceIdT1[b].size();
            m_bndFluxWeightPerTrace[b].assign(numTrace * m_nComp, TData(1.0));
            m_bndGhostIsReflected[b].assign(numTrace, false);
            m_bndCondDiffState[b].assign(numTrace, nullptr);
            m_bndCondDiffStateBlk[b].assign(numTrace, 0u);

            if (!this->HaveBndStorage())
            {
                continue;
            }

            for (size_t t = 0; t < numTrace; ++t)
            {
                unsigned blk;
                size_t offset, compOffset;
                this->m_BndCondOp->GetTraceLocation(this->m_bndTraceIdT1[b][t],
                                                    blk, offset, compOffset);

                // A region no condition claims keeps whatever the session
                // evaluated: an ordinary Dirichlet value, imposed and
                // unweighted, which is what the defaults below give.
                const auto *bndCondOp = this->GetBndCondUpdateOp(blk);

                if (bndCondOp)
                {
                    // Reflective only if something has actually built the
                    // mirrored state in the boundary storage, so the operator
                    // that built it is asked rather than the tag being matched
                    // here. An inviscid Wall or a Symmetry plane becomes
                    // reflective by declaring itself so, with no change to
                    // this file.
                    m_bndGhostIsReflected[b][t] = bndCondOp->GhostIsReflected();

                    // Every plasma energy is an independent component;
                    // preserve each boundary condition's component weights.
                    for (unsigned c = 0; c < m_nComp; ++c)
                        m_bndFluxWeightPerTrace[b][c * numTrace + t] =
                            bndCondOp->BndFluxWeight(blk, c);

                    // A condition acting on the diffusion's exterior state is
                    // called back in the boundary loop of v_Apply(), where
                    // that state exists. Recorded per trace so the loop is a
                    // pointer test rather than a virtual call on every trace
                    // of every block.
                    if (bndCondOp->HasDiffusionExteriorState())
                    {
                        m_bndCondDiffState[b][t] = bndCondOp;
                        m_bndCondDiffStateBlk[b][t] = blk;
                    }
                }
            }
        }
    }

    /// Key this operator is registered under in the OperatorFactory. Defined
    /// by the generated factory boilerplate, not in this header.
    // className - for OperatorFactory
    static std::string className;

    /**
     * @brief Factory creator function.
     *
     * @param expansionList - Expansion list the operator acts on.
     * @param components    - Names of the conserved variables.
     *
     * @return A new operator instance, owned by the caller.
     */
    // instantiation function for CreatorFunction in OperatorFactory
    static std::unique_ptr<MultiRegions::Operator<TData>>
    Instantiate(const MultiRegions::ExpListSharedPtr &expansionList,
                const std::vector<std::string> &components)
    {
        return std::make_unique<FluxTraceAdvDiffOpImpl<
            RiemannKernel, ClosureType, ExecSpace, TData>>(expansionList,
                                                           components);
    }

  protected:
    /// Interleave width of an element-local field, which the gathers and
    /// scatters read and write it at; Field reshapes all its blocks together.
    template <typename TField> NEK_FORCE_INLINE static unsigned Width(TField &f)
    {
        return f.GetBlocks()[0].GetInterleaveWidth();
    }

    /// Spatial dimension of the mesh.
    unsigned m_nDim;
    /// Number of plasma state components, including electron internal energy.
    unsigned m_nComp;

    // workspace
    /// Backing storage for #m_gloDerivT0, #m_gloDerivT1, #m_penFactor and
    /// #m_bndFluxWeight, carved up in the constructor. Anything added to
    /// that carving needs its length added to the size there too.
    LibUtilities::MemoryRegion<TData> m_wspDeriv;
    LibUtilities::MemoryRegion<TData> m_wspGloT1Diff;
    /// Forward-side gradient of the conserved state, packed on the global
    /// trace with `nDim * nComp` components per point.
    TData *m_gloDerivT0;
    /// Backward-side gradient, same layout as #m_gloDerivT0.
    TData *m_gloDerivT1;

    /// Per trace point weight on the energy component of the viscous boundary
    /// flux: zero on an adiabatic wall, one everywhere else. Filled alongside
    /// the penalty factor in the boundary loop.
    TData *m_bndFluxWeight;
    /// Per boundary trace, resolved once from the USERDEFINEDTYPE of the
    /// region it belongs to. Indexed as m_bndTraceIdT1.
    std::vector<std::vector<TData>> m_bndFluxWeightPerTrace;
    /// Per-point geometric penalty factor gathered from the data warehouse.
    TData *m_penFactor;

    /// Whether a boundary trace's exterior state is already a reflection of
    /// the interior about the condition it imposes, as a wall's mirrored
    /// momentum is. Indexed as m_bndTraceIdT1.
    std::vector<std::vector<bool>> m_bndGhostIsReflected;

    /// The condition acting on the diffusion's exterior state at each boundary
    /// trace, or null where none does. Indexed as m_bndTraceIdT1.
    std::vector<std::vector<const SolverCore::BndCondUpdateOp<TData> *>>
        m_bndCondDiffState;

    /// Boundary storage block of each such trace, which the condition needs to
    /// tell its own regions apart. Indexed as m_bndTraceIdT1.
    std::vector<std::vector<unsigned>> m_bndCondDiffStateBlk;

    /**
     * @brief The boundary-only half of a block's gather metadata, resident in
     * the execution space.
     *
     * Where each trace sits in the global trace is common to every trace-flux
     * operator, and the base class mirrors it as
     * TraceFluxOpImpl::m_bndGloTraceDev. What is particular to this operator
     * is the pair below, which follow the boundary conditions rather than the
     * mesh and so are resolved by BuildBndFluxWeights() rather than by the
     * constructor. Hence the separate, lazy upload guarded by
     * #m_bndGatherDataBuilt.
     */
    struct BndGatherDeviceData
    {
        /// Per-component, per-trace diffusion weights (component-major).
        LibUtilities::MemoryRegion<TData> fluxWeight;
        /// Per-trace reflected-ghost flag. The region could hold `bool`; it
        /// holds `unsigned` because it is filled from a `std::vector`, and
        /// `std::vector<bool>` has no contiguous storage to copy from.
        LibUtilities::MemoryRegion<unsigned> ghostIsReflected;
    };

    /// Boundary blocks, indexed as m_bndT0.
    std::vector<BndGatherDeviceData> m_bndGather;
    /// Whether the above has been filled.
    bool m_bndGatherDataBuilt{false};

    /**
     * @brief Make the boundary gather metadata available.
     *
     * Called from DirichletOp(). Pulls BuildBndFluxWeights() forward when it
     * has to: the weights are part of what is uploaded, so they cannot wait
     * for the boundary loop that reads them.
     */
    void EnsureBndGatherData()
    {
        if (m_bndGatherDataBuilt)
        {
            return;
        }

        // An adiabatic wall is a condition on the heat flux, not on the
        // state: with u = 0 at the wall the energy component of the viscous
        // flux is exactly q.n, so suppressing it is the condition. Resolve
        // which traces that applies to once.
        if (m_bndFluxWeightPerTrace.empty())
        {
            BuildBndFluxWeights();
        }

        m_bndGather.clear();
        m_bndGather.resize(this->m_bndT0.size());

        for (unsigned b = 0; b < this->m_bndT0.size(); ++b)
        {
            m_bndGather[b].fluxWeight =
                LibUtilities::MemoryRegion<TData>::template FromVector<
                    MemSpace>(m_bndFluxWeightPerTrace[b]);

            // std::vector<bool> is a bitfield, so it is widened on the host
            // before being copied rather than pointed at.
            std::vector<unsigned> reflected(m_bndGhostIsReflected[b].begin(),
                                            m_bndGhostIsReflected[b].end());
            m_bndGather[b].ghostIsReflected = LibUtilities::MemoryRegion<
                unsigned>::template FromVector<MemSpace>(reflected);
        }

        m_bndGatherDataBuilt = true;
    }

    /**
     * @brief Exterior state the *diffusion* sees on a boundary trace.
     *
     * Zhen Guo's thesis, Ch. 2, "Symmetric interior penalty Galerkin method",
     * Eq. (2.17), writes the interior penalty flux with a single averaging
     * rule, {w} = (w+ + w-)/2, and penalty beta_p = (P+1)^2 / h_T. That rule
     * is applied unchanged on a physical boundary, with the exterior state w-
     * constructed from the boundary condition (Ch. 3.2: "on physical
     * boundaries, there is a complex relation between Q+ and Q-, which are
     * implicitly defined by the physical boundary conditions").
     *
     * A half-and-half average only reproduces the imposed state g if w- is the
     * *reflection* 2g - Q+. A wall ghost already is one - mirroring momentum
     * reflects about zero velocity - but a Dirichlet ghost is g itself, since
     * that same array is what the Riemann solver needs for the advective flux.
     * Pairing g with a half average would impose only half of it and halve the
     * penalty jump with it, so the diffusion is given its own exterior state
     * here: the reflection where the condition is Dirichlet, and m_gloT1
     * unchanged where the ghost already reflects.
     *
     * The same numbers could be reached by sharing one exterior state and
     * varying the average weight per region instead, 1 where the value is
     * imposed and 1/2 where the ghost reflects. Both give {w} = g and
     * jump = 2(g - Q+); the difference is only whether the condition-specific
     * part lives in the ghost or in the weight. Keeping it in the ghost leaves
     * the kernel a direct transcription of Eq. (2.17), with one averaging rule.
     */
    TData *m_gloT1Diff;

    /// Whether the viscous trace flux takes the average-of-normal-fluxes
    /// form, in which the second input of Apply() is each element's outward
    /// normal viscous flux on its traces rather than the gradient.
    bool UseAvgNormalFlux() const
    {
        return this->GetViscousTraceForm() ==
               SolverCore::TraceFluxOp<
                   TData>::ViscousTraceForm::AverageOfNormalFluxes;
    }

    /// The gradient trace the boundary traces read: the second input under
    /// the flux-of-averages form, the field set through
    /// SetBoundaryGradientTrace() otherwise.
    LibUtilities::Field<TData, FieldState::Phys> &
    BoundaryGradientTrace(LibUtilities::Field<TData, FieldState::Phys> &inDeriv)
    {
        if (!UseAvgNormalFlux())
        {
            return inDeriv;
        }
        auto *grad = this->GetBoundaryGradientTrace();
        ASSERTL0(grad || !this->HasBoundaryTraces(),
                 "The average-of-normal-fluxes form needs the gradient trace "
                 "for the boundary traces; see SetBoundaryGradientTrace().");
        return grad ? *grad : inDeriv;
    }

    typename ClosureType::Parameters m_parameters;
    TData m_anomD, m_anomNu, m_anomChi, m_penaltyCoeff;
    TData m_diffusionScale = 1;
    unsigned m_electronEnergyId = std::numeric_limits<unsigned>::max();
    LibUtilities::MemoryRegion<TData> m_wspDiffusion;
    TraceDiffusionWorkspace<TData> m_diffusionWorkspace;

    TraceDiffusionData<ClosureType, TData> DiffusionData()
    {
        ASSERTL0(m_electronEnergyId < m_nComp,
                 "Set the electron-energy component before trace diffusion.");
        return {
            m_parameters,
            m_n_species,
            m_electronEnergyId,
            m_n_ids.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
            m_v_ids.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
            m_e_ids.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
            m_masses.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
            m_charges.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
            m_anomD,
            m_anomNu,
            m_anomChi,
            m_penaltyCoeff,
            m_diffusionScale};
    }

    unsigned int m_n_species = 0;
    LibUtilities::MemoryRegion<unsigned int> m_n_ids;
    LibUtilities::MemoryRegion<unsigned int> m_v_ids;
    LibUtilities::MemoryRegion<unsigned int> m_e_ids;
    LibUtilities::MemoryRegion<TData> m_masses;
    LibUtilities::MemoryRegion<TData> m_charges;

  public:
    void SetElectronEnergyIndex(unsigned id) override
    {
        ASSERTL0(id < m_nComp, "Electron-energy index is outside the state.");
        m_electronEnergyId = id;
    }
    // Must match the scale used to form the input volume normal fluxes.
    void SetDiffusionScale(TData scale) override
    {
        ASSERTL0(std::isfinite(scale), "Diffusion scale must be finite.");
        m_diffusionScale = scale;
    }

    void SetIndices(const std::vector<unsigned int> &n_ids,
                    const std::vector<unsigned int> &v_ids,
                    const std::vector<unsigned int> &e_ids,
                    const std::vector<TData> &masses,
                    const std::vector<TData> &charges) override
    {
        constexpr unsigned invalid = std::numeric_limits<unsigned>::max();
        ASSERTL0(n_ids.size() == v_ids.size() && n_ids.size() == e_ids.size() &&
                     n_ids.size() == masses.size() &&
                     n_ids.size() == charges.size(),
                 "Species metadata vectors must have equal lengths.");
        std::vector<bool> used(m_nComp, false);
        for (size_t spec = 0; spec < n_ids.size(); ++spec)
        {
            ASSERTL0(std::isfinite(masses[spec]) && masses[spec] > 0 &&
                         std::isfinite(charges[spec]) && charges[spec] >= 0,
                     "Invalid species mass or charge.");
            for (unsigned id : {n_ids[spec], v_ids[spec], e_ids[spec]})
                if (id != invalid)
                {
                    ASSERTL0(
                        id < m_nComp && !used[id],
                        "Species component indices must be valid and unique.");
                    used[id] = true;
                }
            ASSERTL0(charges[spec] == 0 ||
                         (n_ids[spec] != invalid && e_ids[spec] != invalid),
                     "Charged species require density and internal energy.");
        }
        m_speciesComponents = std::move(used);
        m_indicesConfigured = true;
        m_n_species = n_ids.size();
        m_n_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(n_ids);
        m_v_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(v_ids);
        m_e_ids = LibUtilities::MemoryRegion<unsigned int>::template FromVector<
            MemSpace, unsigned int>(e_ids);
        m_masses = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(masses);
        m_charges = LibUtilities::MemoryRegion<TData>::template FromVector<
            MemSpace, TData>(charges);
    }

  protected:
    // The base declares each hook twice, for the two- and three-field forms;
    // overriding one form alone would hide the other.
    using SolverCore::TraceFluxOp<TData>::v_Apply;
    using SolverCore::TraceFluxOp<TData>::v_ApplyParallel;
    /**
     * @brief Evaluate the combined inviscid and viscous trace flux on the
     * interior and boundary traces.
     *
     * Dispatches to InteriorOp() and BoundaryOp() for the trace dimension
     * determined at construction. Neither needs anything from another rank,
     * so a caller may run this between BeginParallelExchange() and
     * EndParallelExchange(); the partition traces are v_ApplyParallel()'s,
     * afterwards.
     *
     * @param in      - Physical-space conserved variables.
     * @param inDeriv - Physical-space gradient of @p in, `nDim` components per
     *                  conserved variable; under the average-of-normal-fluxes
     *                  form each element's outward normal viscous flux, one
     *                  component per conserved variable.
     * @param flux    - Trace-space output.
     */
    void v_Apply(LibUtilities::Field<TData, FieldState::Phys> &in,
                 LibUtilities::Field<TData, FieldState::Phys> &inDeriv,
                 LibUtilities::Field<TData, FieldState::Phys> &flux) override
    {
        this->JoinBlockStreams();

        ValidateInputs(in, inDeriv, flux);
        PrepareAuxiliary();
        auto &bndDeriv = BoundaryGradientTrace(inDeriv);
        // The flux is written at the width of the state traces; it is zeroed
        // below, so taking that width needs no reshape.
        flux.SetInterleaveWidth(in);

        // primarily here to ensure all blocks are intialised since
        // the GetPtr() below will only intiialise first block. On the trace
        // stream, which the kernels writing the flux run on.
        if (!this->m_append)
            flux.template Initialize<MemSpace>(0, this->m_opStreamID);

        // The boundary pass calls the conditions' diffusion hook.
        this->SetBndCondStreams();

        // loop over interior and boundary traces and evaluate flux
        // using InteriorOp and BoundaryOp defined below
        switch (this->m_traceDim)
        {
        case 0:
            this->template InteriorOp<0>(in, inDeriv, flux);
            this->template BoundaryOp<0>(in, bndDeriv, flux);
            break;
        case 1:
            this->template InteriorOp<1>(in, inDeriv, flux);
            this->template BoundaryOp<1>(in, bndDeriv, flux);
            break;
        case 2:
            this->template InteriorOp<2>(in, inDeriv, flux);
            this->template BoundaryOp<2>(in, bndDeriv, flux);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Uknown value of m_traceDim");
            break;
        }
        this->ResetBndCondStreams();
        if (this->m_parT0.empty())
            m_auxReady = false;
        this->ReleaseBlockStreams();
    }

    /**
     * @brief Evaluate the combined inviscid and viscous trace flux on the
     * partition-cut traces.
     *
     * Reads the state and gradient the two exchange channels delivered, so it
     * must follow EndParallelExchange(); @p flux must be the field v_Apply()
     * wrote, whose partition slots this fills. Does nothing at all in serial.
     */
    void
    v_ApplyParallel(LibUtilities::Field<TData, FieldState::Phys> &in,
                    LibUtilities::Field<TData, FieldState::Phys> &inDeriv,
                    LibUtilities::Field<TData, FieldState::Phys> &flux) override
    {
        if (this->m_parT0.empty())
        {
            return;
        }

        ValidateInputs(in, inDeriv, flux);
        ASSERTL0(m_auxReady && !m_exchangeActive[0] && !m_exchangeActive[1],
                 "Complete the trace exchanges before ApplyParallel.");
        this->JoinBlockStreams();

        switch (this->m_traceDim)
        {
        case 0:
            this->template ParallelOp<0>(in, inDeriv, flux);
            break;
        case 1:
            this->template ParallelOp<1>(in, inDeriv, flux);
            break;
        case 2:
            this->template ParallelOp<2>(in, inDeriv, flux);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Uknown value of m_traceDim");
            break;
        }
        m_auxReady = false;
        this->ReleaseBlockStreams();
    }

    /**
     * @brief Evaluate the combined flux on every interior trace.
     *
     * For each interior block: gather the conserved state and its gradient on
     * both sides, gather the trace normals and penalty factor, solve the
     * Riemann problem, negate the result and accumulate the viscous flux on
     * top in one pass, and scatter back to both elements.
     *
     * @tparam TRACEDIM Trace dimension; must equal `m_traceDim`.
     * @param  in       Physical-space conserved variables.
     * @param  inDeriv  Physical-space gradient of @p in.
     * @param  flux     Trace-space output.
     */
    template <unsigned TRACEDIM>
    NEK_FORCE_INLINE void
    InteriorOp(LibUtilities::Field<TData, FieldState::Phys> &in,
               LibUtilities::Field<TData, FieldState::Phys> &inDeriv,
               LibUtilities::Field<TData, FieldState::Phys> &flux)
    {
        auto numflux = flux.GetNumComponents();
        const unsigned dim = TRACEDIM + 1;
        const bool avgNormalFlux = UseAvgNormalFlux();
        // Components of the second input: the gradient, or under the
        // average-of-normal-fluxes form the normal viscous flux.
        const unsigned numSecond = avgNormalFlux ? numflux : numflux * m_nDim;
        // The gather scratch is sliced per (trace, component), so it is
        // sized for the larger of the two gathers.
        this->SetWorkSpace(numSecond);

        ASSERTL1(in.GetNumComponents() == numflux,
                 "Assumption that input components are the "
                 "same as the flux components");

        // get field pointers
        auto inPtr = in.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
            this->m_opStreamID);
        auto inDerivPtr =
            inDeriv.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
                this->m_opStreamID);
        auto fluxPtr = flux.GetBlocks()[0].template GetPtr<MemSpace, ReadWrite>(
            this->m_opStreamID);
        auto normPtr = this->m_traceNormals;

        this->template SyncTrailingBlocks<ReadOnly>(in);
        this->template SyncTrailingBlocks<ReadOnly>(inDeriv);
        this->template SyncTrailingBlocks<ReadWrite>(flux);

        // These trace geometry/weight arrays are read-only mesh data, so
        // they live in the data warehouse instead of per-op Field storage.
        // One array spanning every block, one component per block - so it is
        // indexed by the scalar offset below and not by the normals'.
        auto penFactorPtr = this->m_dataWarehouse->template GetData<MemSpace>(
            MultiRegions::IPTraceScalarKey<TData>(
                0, MultiRegions::IPTraceScalarData::IPPenaltyFactor));

        for (unsigned b = 0; b < this->m_intT0.size(); ++b)
        {
            // Step 1: Get global traces associated with block b and put into
            // m_gloT0 and m_gloT1
            this->template SetInteriorParams<TRACEDIM>(b);
            this->template GetInteriorTraces<TRACEDIM>(
                b, numflux, inPtr, this->m_gloT0, this->m_gloT1, Width(in));
            this->template GetInteriorTraces<TRACEDIM>(
                b, numSecond, inDerivPtr, m_gloDerivT0, m_gloDerivT1,
                Width(inDeriv));
            // Step 2: Solve Riemann proble on close packed variables
            // left and right  traces (m_gloT0, m_gloT1),
            auto numBlock = this->m_intT0[b].offset.size();
            auto npTBlock = this->m_intT0[b].m_nTraceXnPtsPad;
            for (unsigned t = 0; t < numBlock; ++t)
            {
                // Single component per block, so not the normals' offset,
                // which is in their coordDim-component layout.
                //
                // Worth checking: the failure this catches - indexing this
                // array with the normals' coordDim-component offset - reads
                // adjacent trace data far more often than anything invalid, so
                // it gives a plausible wrong answer rather than a crash, and
                // only on a mesh with more than one trace block.
                ASSERTL1(
                    this->m_intGloTraceScalarOffset[b][t] + this->m_npTot <=
                        this->m_gloTraceScalarSize,
                    "Single-component trace data is read past the end of "
                    "its block: offset " +
                        std::to_string(this->m_intGloTraceScalarOffset[b][t]) +
                        " + " + std::to_string(this->m_npTot) + " points > " +
                        std::to_string(this->m_gloTraceScalarSize) + ".");
            }

            // gather the normals and the penalty factor
            const auto G =
                this->GetGloTraceOffsetView(this->m_intGloTraceDev[b]);

            SolverCore::detail::GatherGloTraceComponentsKernel<ExecSpace>(
                numBlock, this->m_npTot, npTBlock, dim, normPtr, G,
                this->m_norms, this->m_opStreamID);
            SolverCore::detail::GatherGloTraceScalarKernel<ExecSpace>(
                numBlock, this->m_npTot, penFactorPtr, G, m_penFactor,
                this->m_opStreamID);

            const auto data = DiffusionData();
            auto &aux = m_intAux[b];
            const auto *bTrace =
                aux.b0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID);
            FluxKernelLauncher<RiemannKernel, ExecSpace>(
                npTBlock, m_nDim, m_nComp, data.nSpecies, data.electronEnergyId,
                data.nIds, data.vIds, data.eIds, data.masses, data.charges,
                this->m_norms, bTrace,
                aux.j0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                aux.j1.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                this->m_gloT0, this->m_gloT1, this->m_flux, m_advectionScale,
                this->m_opStreamID);

            if (avgNormalFlux)
            {
                // The second input holds each side's outward normal viscous
                // flux; the kernel averages them and adds the penalty.
                AvgNormalFluxTraceKernel<ExecSpace, ClosureType>(
                    DiffusionData(), npTBlock, m_nDim, m_nComp, this->m_norms,
                    m_penFactor, bTrace, this->m_gloT0, this->m_gloT1,
                    m_gloDerivT0, m_gloDerivT1, this->m_flux,
                    m_diffusionWorkspace, this->m_opStreamID);
            }
            else
            {
                DiffuseTraceFluxKernel<ExecSpace, ClosureType, true>(
                    DiffusionData(), npTBlock, m_nDim, m_nComp, this->m_norms,
                    m_penFactor, static_cast<const TData *>(nullptr), bTrace,
                    this->m_gloT0, this->m_gloT1, m_gloDerivT0, m_gloDerivT1,
                    this->m_flux, m_diffusionWorkspace, this->m_opStreamID);
            }
            // Step 3: Replace flux values from (m_gloT0,m_gloT1) into
            // the local flux storage
            if (this->m_append)
            {
                this->template InterpBackInteriorFlux<TRACEDIM, true>(
                    b, numflux, fluxPtr, Width(flux));
            }
            else
            {
                this->template InterpBackInteriorFlux<TRACEDIM, false>(
                    b, numflux, fluxPtr, Width(flux));
            }
        }
    }

    /**
     * @brief Evaluate the combined flux on every trace cut by the partitioner.
     *
     * InteriorOp() for traces whose backward element is on another rank, and
     * identical to it in every respect the physics can see. A parallel trace
     * is an *interior* trace: its exterior state and gradient are a real
     * neighbouring element's, delivered on exchange channels 0 and 1, so the
     * Riemann solve is the interior one and DiffuseTraceFluxKernel() is
     * instantiated with `IsInterior = true` - average weights of one half,
     * jump weights of one, no boundary energy-flux weight. Routing these
     * through the boundary form would reflect a perfectly good neighbour
     * state about itself and return something finite, smooth and wrong.
     *
     * The result is scattered to the local element alone; the neighbour
     * computes the same flux from the same four traces and scatters it to
     * the other element itself.
     *
     * Returns immediately when there are no such traces, which is every case
     * in serial.
     *
     * @tparam TRACEDIM Trace dimension; must equal `m_traceDim`.
     * @param  in       Physical-space conserved variables.
     * @param  inDeriv  Physical-space gradient of @p in.
     * @param  flux     Trace-space output.
     */
    template <unsigned TRACEDIM>
    NEK_FORCE_INLINE void
    ParallelOp(LibUtilities::Field<TData, FieldState::Phys> &in,
               LibUtilities::Field<TData, FieldState::Phys> &inDeriv,
               LibUtilities::Field<TData, FieldState::Phys> &flux)
    {
        if (this->m_parT0.empty())
        {
            return;
        }

        auto numflux = flux.GetNumComponents();
        const unsigned dim = TRACEDIM + 1;
        const bool avgNormalFlux = UseAvgNormalFlux();
        // Components of the second input: the gradient, or under the
        // average-of-normal-fluxes form the normal viscous flux.
        const unsigned numSecond = avgNormalFlux ? numflux : numflux * m_nDim;

        // The gather scratch is sliced per (trace, component), so it is
        // sized for the larger of the two gathers.
        this->SetWorkSpace(numSecond);

        ASSERTL1(in.GetNumComponents() == numflux,
                 "Assumption that input components are the "
                 "same as the flux components");

        // get field pointers
        auto inPtr = in.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
            this->m_opStreamID);
        auto inDerivPtr =
            inDeriv.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
                this->m_opStreamID);
        auto fluxPtr = flux.GetBlocks()[0].template GetPtr<MemSpace, ReadWrite>(
            this->m_opStreamID);
        auto normPtr = this->m_traceNormals;

        this->template SyncTrailingBlocks<ReadOnly>(in);
        this->template SyncTrailingBlocks<ReadOnly>(inDeriv);
        this->template SyncTrailingBlocks<ReadWrite>(flux);

        auto penFactorPtr = this->m_dataWarehouse->template GetData<MemSpace>(
            MultiRegions::IPTraceScalarKey<TData>(
                0, MultiRegions::IPTraceScalarData::IPPenaltyFactor));

        for (unsigned b = 0; b < this->m_parT0.size(); ++b)
        {
            // Step 1: Get global traces associated with block b and put into
            // m_gloT0 and m_gloT1; the backward sides come from the receive
            // buffers of channels 0 and 1.
            this->template SetParallelParams<TRACEDIM>(b);
            this->template GetParallelTraces<TRACEDIM>(
                b, numflux, inPtr, this->m_gloT0, this->m_gloT1, 0, Width(in));
            this->template GetParallelTraces<TRACEDIM>(
                b, numSecond, inDerivPtr, m_gloDerivT0, m_gloDerivT1, 1,
                Width(inDeriv));

            // Step 2: gather the normals and penalty factor for these traces
            auto numBlock = this->m_parT0[b].offset.size();
            auto npTBlock = this->m_parT0[b].m_nTraceXnPtsPad;
            for (unsigned t = 0; t < numBlock; ++t)
            {
                // Single component per block, so not the normals' offset; see
                // the note on the same read in InteriorOp().
                ASSERTL1(
                    this->m_parGloTraceScalarOffset[b][t] + this->m_npTot <=
                        this->m_gloTraceScalarSize,
                    "Single-component trace data is read past the end of "
                    "its block: offset " +
                        std::to_string(this->m_parGloTraceScalarOffset[b][t]) +
                        " + " + std::to_string(this->m_npTot) + " points > " +
                        std::to_string(this->m_gloTraceScalarSize) + ".");
            }

            const auto G =
                this->GetGloTraceOffsetView(this->m_parGloTraceDev[b]);

            SolverCore::detail::GatherGloTraceComponentsKernel<ExecSpace>(
                numBlock, this->m_npTot, npTBlock, dim, normPtr, G,
                this->m_norms, this->m_opStreamID);
            SolverCore::detail::GatherGloTraceScalarKernel<ExecSpace>(
                numBlock, this->m_npTot, penFactorPtr, G, m_penFactor,
                this->m_opStreamID);

            const auto data = DiffusionData();
            auto &aux = m_parAux[b];
            const auto *bTrace =
                aux.b0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID);
            FluxKernelLauncher<RiemannKernel, ExecSpace>(
                npTBlock, m_nDim, m_nComp, data.nSpecies, data.electronEnergyId,
                data.nIds, data.vIds, data.eIds, data.masses, data.charges,
                this->m_norms, bTrace,
                aux.j0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                aux.j1.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                this->m_gloT0, this->m_gloT1, this->m_flux, m_advectionScale,
                this->m_opStreamID);

            if (avgNormalFlux)
            {
                // The second input holds each side's outward normal viscous
                // flux; the kernel averages them and adds the penalty.
                AvgNormalFluxTraceKernel<ExecSpace, ClosureType>(
                    DiffusionData(), npTBlock, m_nDim, m_nComp, this->m_norms,
                    m_penFactor, bTrace, this->m_gloT0, this->m_gloT1,
                    m_gloDerivT0, m_gloDerivT1, this->m_flux,
                    m_diffusionWorkspace, this->m_opStreamID);
            }
            else
            {
                DiffuseTraceFluxKernel<ExecSpace, ClosureType, true>(
                    DiffusionData(), npTBlock, m_nDim, m_nComp, this->m_norms,
                    m_penFactor, static_cast<const TData *>(nullptr), bTrace,
                    this->m_gloT0, this->m_gloT1, m_gloDerivT0, m_gloDerivT1,
                    this->m_flux, m_diffusionWorkspace, this->m_opStreamID);
            }

            // Step 3: Replace flux values into the local flux storage, this
            // side only.
            if (this->m_append)
            {
                this->template InterpBackParallelFlux<TRACEDIM, true>(
                    b, numflux, fluxPtr, Width(flux));
            }
            else
            {
                this->template InterpBackParallelFlux<TRACEDIM, false>(
                    b, numflux, fluxPtr, Width(flux));
            }
        }
    }

    /**
     * @brief Evaluate the combined flux on every Dirichlet boundary trace.
     *
     * As InteriorOp(), with the exterior conserved state taken from the
     * boundary condition operator and the exterior gradient taken equal to the
     * interior one, since the boundary condition operator holds no gradient
     * data. DiffuseTraceFluxKernel() is instantiated in its boundary form.
     *
     * @tparam TRACEDIM Trace dimension; must equal `m_traceDim`.
     * @param  in       Physical-space conserved variables.
     * @param  inDeriv  Physical-space gradient of @p in.
     * @param  flux     Trace-space output.
     */
    template <unsigned TRACEDIM>
    NEK_FORCE_INLINE void
    BoundaryOp(LibUtilities::Field<TData, FieldState::Phys> &in,
               LibUtilities::Field<TData, FieldState::Phys> &inDeriv,
               LibUtilities::Field<TData, FieldState::Phys> &flux)
    {
        // Nothing to read or write without a boundary block; the gradient
        // trace in particular may never have been formed on this rank.
        if (this->m_bndT0.empty())
        {
            return;
        }

        auto numflux = flux.GetNumComponents();
        const unsigned dim = TRACEDIM + 1;

        // The gather scratch is sliced per (trace, component), so it is
        // sized for the gradient gather, the boundary traces reading the
        // gradient under either viscous trace form.
        this->SetWorkSpace(numflux * m_nDim);

        ASSERTL1(in.GetNumComponents() == numflux,
                 "Assumption that input components are the "
                 "same as the flux components");

        // get field pointers
        auto inPtr = in.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
            this->m_opStreamID);
        auto inDerivPtr =
            inDeriv.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
                this->m_opStreamID);
        auto fluxPtr = flux.GetBlocks()[0].template GetPtr<MemSpace, ReadWrite>(
            this->m_opStreamID);
        auto normPtr = this->m_traceNormals;

        this->template SyncTrailingBlocks<ReadOnly>(in);
        this->template SyncTrailingBlocks<ReadOnly>(inDeriv);
        this->template SyncTrailingBlocks<ReadWrite>(flux);

        // These trace geometry/weight arrays are read-only mesh data, so
        // they live in the data warehouse instead of per-op Field storage.
        // These trace geometry/weight arrays are read-only mesh data, so
        // they live in the data warehouse instead of per-op Field storage.
        // One array spanning every block, one component per block - so it is
        // indexed by the scalar offset below and not by the normals'.
        auto penFactorPtr = this->m_dataWarehouse->template GetData<MemSpace>(
            MultiRegions::IPTraceScalarKey<TData>(
                0, MultiRegions::IPTraceScalarData::IPPenaltyFactor));

        EnsureBndGatherData();

        for (unsigned b = 0; b < this->m_bndT0.size(); ++b)
        {
            // Step 1: Get global trace 't' and put into m_gloT0 and
            // m_gloT1
            this->template SetBoundaryParams<TRACEDIM>(b);
            this->template GetBoundaryTraces<TRACEDIM>(
                b, numflux, inPtr, this->m_gloT0, this->m_gloT1, Width(in));

            auto numBlock = this->m_bndT0[b].offset.size();
            auto npTBlock = this->m_bndT0[b].m_nTraceXnPtsPad;

            // The normals are gathered before the gradient traces: a Neumann
            // condition - the adiabatic wall's dE/dn - needs them to resolve
            // the normal component of the gradient.
            for (unsigned t = 0; t < numBlock; ++t)
            {
                // Single component per block, so not the normals' offset,
                // which is in their coordDim-component layout.
                //
                // Worth checking: the failure this catches - indexing this
                // array with the normals' coordDim-component offset - reads
                // adjacent trace data far more often than anything invalid, so
                // it gives a plausible wrong answer rather than a crash, and
                // only on a mesh with more than one trace block.
                ASSERTL1(
                    this->m_bndGloTraceScalarOffset[b][t] + this->m_npTot <=
                        this->m_gloTraceScalarSize,
                    "Single-component trace data is read past the end of "
                    "its block: offset " +
                        std::to_string(this->m_bndGloTraceScalarOffset[b][t]) +
                        " + " + std::to_string(this->m_npTot) + " points > " +
                        std::to_string(this->m_gloTraceScalarSize) + ".");
            }

            // gather the normals, the penalty factor, and the energy weight,
            // which is one value per trace rather than one per point
            const auto G =
                this->GetGloTraceOffsetView(this->m_bndGloTraceDev[b]);

            SolverCore::detail::GatherGloTraceComponentsKernel<ExecSpace>(
                numBlock, this->m_npTot, npTBlock, dim, normPtr, G,
                this->m_norms, this->m_opStreamID);
            SolverCore::detail::GatherGloTraceScalarKernel<ExecSpace>(
                numBlock, this->m_npTot, penFactorPtr, G, m_penFactor,
                this->m_opStreamID);
            for (unsigned c = 0; c < m_nComp; ++c)
                SolverCore::detail::BroadcastGloTraceValueKernel<ExecSpace>(
                    numBlock, this->m_npTot,
                    m_bndGather[b]
                            .fluxWeight.template GetPtr<MemSpace, ReadOnly>(
                                this->m_opStreamID) +
                        c * numBlock,
                    m_bndFluxWeight + c * npTBlock, this->m_opStreamID);

            // Build the exterior state the diffusion sees, which is not always
            // the one the Riemann solver sees.
            //
            // Every average across a trace is the half-and-half {w} =
            // (w+ + w-)/2 of Zhen-Guo Yan, "Efficient implicit spectral/hp
            // element DG techniques for compressible flows", PhD thesis,
            // Imperial College London, 2021, Ch. 2, Eq. (2.17), on boundaries
            // as well as in the interior. That is only correct if the exterior
            // state is a ghost
            // whose average with the interior is the value we mean to impose,
            // and boundary regions come in two kinds:
            //
            //   reflective  - a wall or symmetry plane, where BndCondWallCFEOp
            //                 has already built the ghost by negating momentum
            //                 in the boundary storage. m_gloT1 is (rho+,
            //                 -rhou+, E+), the average is the no-slip state,
            //                 and the ghost stands as it is for both operators.
            //
            //   imposed     - an ordinary session Dirichlet region, where
            //                 m_gloT1 holds the target value g. The Riemann
            //                 solver wants g, so m_gloT1 must keep it; the
            //                 diffusion wants the viscous flux evaluated at g,
            //                 which under a half-and-half average means the
            //                 ghost 2g - Q+. Hence the separate array.
            //
            // The same distinction could be carried as a per-region average
            // weight, 1/2 on reflective regions and 1 elsewhere, applied to the
            // unmodified boundary value: reflecting the ghost and averaging is
            // the same as not reflecting it and taking the boundary value
            // alone. Carrying it in the state keeps one averaging rule in the
            // kernels, so a new boundary condition declares which kind of
            // exterior state it builds instead of threading a weight through
            // the flux evaluation.
            //
            // Seed the whole padded block, not only the points in use: the
            // kernel walks every point of it.
            Nektar::Math::copyKernel<ExecSpace>(numflux * npTBlock,
                                                this->m_gloT1, m_gloT1Diff,
                                                this->m_opStreamID);

            ReflectDiffusionGhostStateKernel<ExecSpace>(
                numBlock, this->m_npTot, npTBlock, numflux,
                m_bndGather[b]
                    .ghostIsReflected.template GetPtr<MemSpace, ReadOnly>(
                        this->m_opStreamID),
                this->m_gloT0, this->m_gloT1, m_gloT1Diff, this->m_opStreamID);

            // Any condition that acts on the diffusion's exterior state
            // and on nothing else - an isothermal wall's temperature, which
            // belongs in the viscous flux and must not reach the Riemann
            // solver. The condition owns the physics; this
            // operator owns only the layout, so it hands over a window onto
            // one trace and the stride between components.
            for (unsigned t = 0; t < numBlock; ++t)
            {
                auto *bndCondOp = m_bndCondDiffState[b][t];

                if (!bndCondOp)
                {
                    continue;
                }

                const size_t o = this->m_npTot * t;

                bndCondOp->ApplyDiffusionExteriorState(
                    m_bndCondDiffStateBlk[b][t], m_gloT1Diff + o,
                    this->m_gloT0 + o, this->m_npTot, npTBlock);
            }

            // The gradient carries any Neumann conditions - for a compressible
            // wall that is dE/dn, which with u = 0 at the wall is the
            // adiabatic dT/dn. Components with any other condition keep the
            // exterior gradient equal to the interior, as before.
            this->template GetNeumannBoundaryTraces<TRACEDIM>(
                b, numflux, m_nDim, inDerivPtr, this->m_norms, m_gloDerivT0,
                m_gloDerivT1, Width(inDeriv));

            // Step 2: Solve Riemann proble on  close packed variables
            // left and right  traces (m_gloT0, m_gloT1),
            const auto data = DiffusionData();
            auto &aux = m_bndAux[b];
            const auto *bTrace =
                aux.b0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID);
            FluxKernelLauncher<RiemannKernel, ExecSpace>(
                npTBlock, m_nDim, m_nComp, data.nSpecies, data.electronEnergyId,
                data.nIds, data.vIds, data.eIds, data.masses, data.charges,
                this->m_norms, bTrace,
                aux.j0.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                aux.j1.template GetPtr<MemSpace, ReadOnly>(this->m_opStreamID),
                this->m_gloT0, this->m_gloT1, this->m_flux, m_advectionScale,
                this->m_opStreamID);

            DiffuseTraceFluxKernel<ExecSpace, ClosureType, false>(
                DiffusionData(), npTBlock, m_nDim, m_nComp, this->m_norms,
                m_penFactor, m_bndFluxWeight, bTrace, this->m_gloT0,
                m_gloT1Diff, m_gloDerivT0, m_gloDerivT1, this->m_flux,
                m_diffusionWorkspace, this->m_opStreamID);

            // Step 3: Replace flux values from (m_gloT0,m_gloT1) into
            // the local flux storage
            if (this->m_append)
            {
                this->template InterpBackDirichletFlux<TRACEDIM, true>(
                    b, numflux, fluxPtr, Width(flux));
            }
            else
            {
                this->template InterpBackDirichletFlux<TRACEDIM, false>(
                    b, numflux, fluxPtr, Width(flux));
            }
        }
    }

  public:
    void
    SetMagneticTrace(LibUtilities::Field<TData, FieldState::Phys> &b) override
    {
        ASSERTL0(!m_exchangeActive[0] && !m_exchangeActive[1],
                 "Cannot replace magnetic traces during an exchange.");
        ASSERTL0(b.GetNumComponents() == 3,
                 "The magnetic unit-vector trace needs three components.");
        m_magneticTrace = &b;
        m_auxReady = false;
    }
    void SetParallelCurrentTrace(
        LibUtilities::Field<TData, FieldState::Phys> &current) override
    {
        ASSERTL0(!m_exchangeActive[0] && !m_exchangeActive[1],
                 "Cannot replace current traces during an exchange.");
        ASSERTL0(current.GetNumComponents() == 1,
                 "The parallel number-current trace needs one component.");
        m_currentTrace = &current;
        m_auxReady = false;
    }
    void SetAdvectionScale(TData scale) override
    {
        ASSERTL0(std::isfinite(scale), "Advection scale must be finite.");
        m_advectionScale = scale;
    }

  protected:
    struct AuxiliaryBlock
    {
        LibUtilities::MemoryRegion<TData> b0, b1, j0, j1;
    };
    std::vector<AuxiliaryBlock> m_intAux, m_parAux, m_bndAux;
    LibUtilities::Field<TData, FieldState::Phys> *m_magneticTrace = nullptr;
    LibUtilities::Field<TData, FieldState::Phys> *m_currentTrace = nullptr;
    bool m_auxReady = false;
    bool m_exchangeActive[2] = {false, false};
    bool m_indicesConfigured = false;
    std::vector<bool> m_speciesComponents;
    TData m_advectionScale = 1;

    template <typename Blocks>
    void AllocateAuxiliary(const Blocks &blocks,
                           std::vector<AuxiliaryBlock> &aux)
    {
        aux.resize(blocks.size());
        for (size_t b = 0; b < blocks.size(); ++b)
        {
            const size_t npts = blocks[b].m_nTraceXnPtsPad;
            aux[b].b0 = LibUtilities::MemoryRegion<TData>(3 * npts);
            aux[b].b1 = LibUtilities::MemoryRegion<TData>(3 * npts);
            aux[b].j0 = LibUtilities::MemoryRegion<TData>(npts);
            aux[b].j1 = LibUtilities::MemoryRegion<TData>(npts);
            for (auto *region :
                 {&aux[b].b0, &aux[b].b1, &aux[b].j0, &aux[b].j1})
                region->template Initialize<MemSpace>(0);
        }
    }

    void ValidateInputs(LibUtilities::Field<TData, FieldState::Phys> &in,
                        LibUtilities::Field<TData, FieldState::Phys> &second,
                        LibUtilities::Field<TData, FieldState::Phys> &flux)
    {
        ASSERTL0(
            m_indicesConfigured && m_magneticTrace && m_currentTrace &&
                m_electronEnergyId < m_nComp &&
                !m_speciesComponents[m_electronEnergyId],
            "Configure species, electron energy, magnetic and current traces.");
        const unsigned numSecond =
            UseAvgNormalFlux() ? m_nComp : m_nComp * m_nDim;
        ASSERTL0(in.GetNumComponents() == m_nComp &&
                     second.GetNumComponents() == numSecond &&
                     flux.GetNumComponents() == m_nComp && &in != &flux &&
                     &second != &flux && m_magneticTrace != &flux &&
                     m_currentTrace != &flux,
                 "Trace state, gradient/normal flux, or output is invalid.");
        auto &boundaryGradient = BoundaryGradientTrace(second);
        if (this->HasBoundaryTraces())
            ASSERTL0(boundaryGradient.GetNumComponents() == m_nComp * m_nDim,
                     "Boundary diffusion requires a full gradient trace.");
        for (auto *field : {&second, &flux, m_magneticTrace, m_currentTrace,
                            &boundaryGradient})
        {
            ASSERTL0(field->GetBlocks().size() == in.GetBlocks().size() &&
                         field->GetNumHomoModes() == in.GetNumHomoModes(),
                     "Element-local trace field block layouts differ.");
            for (size_t b = 0; b < in.GetBlocks().size(); ++b)
            {
                const auto &reference = in.GetBlocks()[b];
                const auto &block = field->GetBlocks()[b];
                ASSERTL0(block.CompSize() == reference.CompSize() &&
                             block.GetNumElements() ==
                                 reference.GetNumElements(),
                         "Element-local trace block sizes differ.");
            }
        }
        // Changing layout for append must preserve the existing output values.
        // The Field reshape handles this before any pointer is obtained.
    }

    template <unsigned TRACEDIM>
    void GatherAuxiliary(LibUtilities::Field<TData, FieldState::Phys> &field,
                         const bool magnetic)
    {
        const unsigned numComp = field.GetNumComponents();
        this->SetWorkSpace(numComp);
        this->JoinBlockStreams();
        auto ptr = field.GetBlocks()[0].template GetPtr<MemSpace, ReadOnly>(
            this->m_opStreamID);
        this->template SyncTrailingBlocks<ReadOnly>(field);
        const unsigned width = Width(field);
        const auto pointers = [&](AuxiliaryBlock &aux)
        {
            return std::make_pair(
                (magnetic ? aux.b0 : aux.j0)
                    .template GetPtr<MemSpace, ReadWrite>(this->m_opStreamID),
                (magnetic ? aux.b1 : aux.j1)
                    .template GetPtr<MemSpace, ReadWrite>(this->m_opStreamID));
        };
        for (unsigned b = 0; b < m_intAux.size(); ++b)
        {
            this->template SetInteriorParams<TRACEDIM>(b);
            const auto [fwd, bwd] = pointers(m_intAux[b]);
            this->template GetInteriorTraces<TRACEDIM>(b, numComp, ptr, fwd,
                                                       bwd, width);
        }
        for (unsigned b = 0; b < m_bndAux.size(); ++b)
        {
            this->template SetBoundaryParams<TRACEDIM>(b);
            const auto [fwd, bwd] = pointers(m_bndAux[b]);
            // Auxiliary fields have no entries in the plasma-state BC storage.
            // Continue the prescribed current and magnetic field to the ghost.
            this->template GetBoundaryTraces<TRACEDIM, false>(b, numComp, ptr,
                                                              fwd, bwd, width);
        }
        if (!this->m_parT0.empty())
        {
            // Use the base's width-aware pack, then consume its receive buffer
            // before channel 0 is reused for the next auxiliary field/state.
            SolverCore::detail::TraceFluxOpImpl<
                ExecSpace, TData,
                FluxTraceAdvDiffOp<TData>>::v_BeginParallelExchange(field, 0);
            SolverCore::detail::TraceFluxOpImpl<
                ExecSpace, TData,
                FluxTraceAdvDiffOp<TData>>::v_EndParallelExchange();
            for (unsigned b = 0; b < m_parAux.size(); ++b)
            {
                this->template SetParallelParams<TRACEDIM>(b);
                const auto [fwd, bwd] = pointers(m_parAux[b]);
                this->template GetParallelTraces<TRACEDIM>(b, numComp, ptr, fwd,
                                                           bwd, 0, width);
            }
        }
    }

    void PrepareAuxiliary()
    {
        if (m_auxReady)
            return;
        ASSERTL0(m_magneticTrace && m_currentTrace && !m_exchangeActive[0] &&
                     !m_exchangeActive[1],
                 "Auxiliary traces must be gathered before state/gradient "
                 "exchanges.");
        switch (this->m_traceDim)
        {
        case 0:
            GatherAuxiliary<0>(*m_magneticTrace, true);
            GatherAuxiliary<0>(*m_currentTrace, false);
            break;
        case 1:
            GatherAuxiliary<1>(*m_magneticTrace, true);
            GatherAuxiliary<1>(*m_currentTrace, false);
            break;
        case 2:
            GatherAuxiliary<2>(*m_magneticTrace, true);
            GatherAuxiliary<2>(*m_currentTrace, false);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "Unsupported trace dimension.");
        }
        m_auxReady = true;
    }

    void
    v_BeginParallelExchange(LibUtilities::Field<TData, FieldState::Phys> &trace,
                            const unsigned chan) override
    {
        ASSERTL0(chan < 2 && !m_exchangeActive[chan],
                 "State and gradient use exchange channels 0 and 1.");
        PrepareAuxiliary();
        SolverCore::detail::TraceFluxOpImpl<
            ExecSpace, TData,
            FluxTraceAdvDiffOp<TData>>::v_BeginParallelExchange(trace, chan);
        if (!this->m_parT0.empty())
            m_exchangeActive[chan] = true;
    }
    void v_EndParallelExchange() override
    {
        SolverCore::detail::TraceFluxOpImpl<
            ExecSpace, TData,
            FluxTraceAdvDiffOp<TData>>::v_EndParallelExchange();
        m_exchangeActive[0] = m_exchangeActive[1] = false;
    }
};

} // namespace PENKNIFE::detail
