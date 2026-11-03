#ifndef MAGNETIC_FIELD_OP_IMPL
#define MAGNETIC_FIELD_OP_IMPL

#include "MagneticFieldKernels.hpp"
#include "MagneticFieldOp.hpp"
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <MultiRegions/ElmtOps/BwdTrans/BwdTransOp.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>

namespace PENKNIFE::detail
{

template <typename ExecSpace, typename TData>
class MagneticFieldOpImpl : public MagneticFieldOp<TData>
{
    using MemSpace = typename ExecSpace::memory_space;
    using Points   = Array<OneD, Array<OneD, NekDouble>>;

public:
    MagneticFieldOpImpl(const MultiRegions::ExpListSharedPtr &expansionList,
                        const std::vector<std::string> &components)
        : MagneticFieldOp<TData>(expansionList, components),
          m_session(expansionList->GetSession()), m_expansion(expansionList),
          m_ptsIO(m_session->GetComm())
    {
        m_dimension = expansionList->GetCoordim(0);
        m_npoints   = expansionList->GetTotPoints();
        m_coords    = Points(3);
        m_values    = Points(3);
        for (unsigned d = 0; d < 3; ++d)
        {
            m_coords[d] = Array<OneD, NekDouble>(m_npoints, 0.0);
            m_values[d] = Array<OneD, NekDouble>(m_npoints, 0.0);
        }
        expansionList->GetCoords(m_coords[0], m_coords[1], m_coords[2]);

        if (m_session->DefinesFunction("MagneticMeanField"))
        {
            m_functionName = "MagneticMeanField";
            if (m_dimension == 3)
            {
                const auto kind =
                    m_session->GetFunctionType(m_functionName, "Bx");
                this->type = kind == LibUtilities::eFunctionTypeFile
                                 ? MagneticFieldOp<TData>::field_type::mean_file
                                 : MagneticFieldOp<TData>::field_type::mean_exp;
                for (const auto &name : m_names)
                {
                    ASSERTL0(m_session->GetFunctionType(m_functionName, name) ==
                                 kind,
                             "A 3D MagneticMeanField must use either three "
                             "expressions or one axisymmetric point file.");
                    if (kind == LibUtilities::eFunctionTypeFile)
                        ASSERTL0(m_session->GetFunctionFilename(m_functionName,
                                                                name) ==
                                     m_session->GetFunctionFilename(
                                         m_functionName, "Bx"),
                                 "Mean-field components must use the same "
                                 "point file.");
                }
                if (this->type == MagneticFieldOp<TData>::field_type::mean_exp)
                {
                    m_radius = Array<OneD, NekDouble>(m_npoints);
                    m_phi    = Array<OneD, NekDouble>(m_npoints);
                    for (size_t p = 0; p < m_npoints; ++p)
                    {
                        m_radius[p] =
                            std::hypot(m_coords[0][p], m_coords[2][p]);
                        m_phi[p] = std::atan2(m_coords[2][p], m_coords[0][p]);
                    }
                }
            }
        }
        else
        {
            ASSERTL0(m_session->DefinesFunction("MagneticField"),
                     "Define MagneticField or MagneticMeanField. "
                     "The CurrentDensity solve will be migrated separately.");
            m_functionName = "MagneticField";
        }
    }

    static std::string className;

    static std::unique_ptr<MultiRegions::Operator<TData>> Instantiate(
        const MultiRegions::ExpListSharedPtr &expansionList,
        const std::vector<std::string> &components)
    {
        return std::make_unique<MagneticFieldOpImpl<ExecSpace, TData>>(
            expansionList, components);
    }

protected:
    void v_Apply(TData time) override
    {
        Read(time);
        // Only file/expression evaluation and point interpolation use host
        // arrays. The solution and derived values live in Fields.
        std::vector<TData> values(3 * m_npoints);
        for (unsigned d = 0; d < 3; ++d)
            std::copy_n(m_values[d].data(), m_npoints,
                        values.data() + d * m_npoints);
        this->m_B->template Initialize<MemSpace>(TData(0));
        this->m_B->template CopyVector<MemSpace>(values);
        this->m_b_unit->SetInterleaveWidth(*this->m_B);
        this->m_mag_B->SetInterleaveWidth(*this->m_B);
        this->m_b_unit->template Initialize<MemSpace>(TData(0));
        this->m_mag_B->template Initialize<MemSpace>(TData(0));

        ASSERTL0(this->m_B->GetBlocks().size() ==
                         this->m_b_unit->GetBlocks().size() &&
                     this->m_B->GetBlocks().size() ==
                         this->m_mag_B->GetBlocks().size(),
                 "Magnetic field outputs must share a block layout.");
        for (unsigned blk = 0; blk < this->m_B->GetBlocks().size(); ++blk)
        {
            const unsigned streamID = blk + 1;
            auto &Bblock            = this->m_B->GetBlocks()[blk];
            auto &unit              = this->m_b_unit->GetBlocks()[blk];
            auto &mag               = this->m_mag_B->GetBlocks()[blk];
            ASSERTL0(Bblock.CompSize() == unit.CompSize() &&
                         Bblock.CompSize() == mag.CompSize(),
                     "Magnetic field outputs must share a block layout.");
            const size_t stride =
                Bblock.CompSize() * this->m_B->GetNumHomoModes();
            MagneticFieldKernel<ExecSpace>(
                stride, stride, stride,
                Bblock.template GetPtr<MemSpace, ReadOnly>(streamID),
                unit.template GetPtr<MemSpace, WriteOnly>(streamID),
                mag.template GetPtr<MemSpace, WriteOnly>(streamID), streamID);
        }
    }

private:
    static inline const std::array<std::string, 3> m_names = {"Bx", "By", "Bz"};
    unsigned m_dimension;
    size_t m_npoints;
    LibUtilities::SessionReaderSharedPtr m_session;
    MultiRegions::ExpListSharedPtr m_expansion;
    LibUtilities::PtsIO m_ptsIO;
    std::string m_functionName;
    Points m_coords, m_values;
    Array<OneD, NekDouble> m_radius, m_phi;
    LibUtilities::PtsFieldSharedPtr m_meanInput, m_meanOutput;
    LibUtilities::Interpolator m_meanInterp{LibUtilities::eShepard};
    Points m_meanCoordinates;
    static constexpr unsigned m_increments = 64;

    void Read(TData time)
    {
        if (this->type == MagneticFieldOp<TData>::field_type::mean_file)
        {
            ReadMeanFile();
            return;
        }
        if (this->type == MagneticFieldOp<TData>::field_type::mean_exp)
        {
            // The expression variables x,y,z denote R,Z,phi in this mode.
            for (unsigned d = 0; d < 3; ++d)
                m_session->GetFunction(m_functionName, m_names[d])
                    ->Evaluate(m_radius, m_coords[1], m_phi, time, m_values[d]);
            for (size_t p = 0; p < m_npoints; ++p)
            {
                const double Br = m_values[0][p], Bphi = m_values[2][p];
                m_values[0][p] =
                    std::cos(m_phi[p]) * Br - std::sin(m_phi[p]) * Bphi;
                m_values[2][p] =
                    std::sin(m_phi[p]) * Br + std::cos(m_phi[p]) * Bphi;
            }
            return;
        }

        // Ordinary Cartesian functions may mix expressions and files.
        std::map<std::string, LibUtilities::PtsFieldSharedPtr> pointFiles;
        for (unsigned d = 0; d < 3; ++d)
        {
            const auto &name = m_names[d];
            if (m_session->GetFunctionType(m_functionName, name) ==
                LibUtilities::eFunctionTypeExpression)
            {
                m_session->GetFunction(m_functionName, name)
                    ->Evaluate(m_coords[0], m_coords[1], m_coords[2], time,
                               m_values[d]);
                continue;
            }
            const auto filename =
                m_session->GetFunctionFilename(m_functionName, name);
            auto variable =
                m_session->GetFunctionFilenameVariable(m_functionName, name);
            if (variable.empty())
                variable = name;
            const auto dot = filename.find_last_of('.');
            if (dot != std::string::npos && filename.substr(dot) == ".pts")
            {
                auto &input = pointFiles[filename];
                if (!input)
                    m_ptsIO.Import(filename, input);
                auto names = input->GetFieldNames();
                const auto field =
                    std::find(names.begin(), names.end(), variable);
                ASSERTL0(field != names.end(),
                         "Missing magnetic point-file field " + variable);
                Points points;
                input->GetPts(points);
                Points target(input->GetDim() + names.size());
                for (unsigned j = 0; j < input->GetDim(); ++j)
                    target[j] = m_coords[j];
                for (unsigned j = 0; j < names.size(); ++j)
                    target[input->GetDim() + j] =
                        Array<OneD, NekDouble>(m_npoints);
                auto output =
                    MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                        input->GetDim(), names, target);
                LibUtilities::Interpolator interp(LibUtilities::eShepard);
                interp.CalcWeights(input, output);
                interp.Interpolate(input, output);
                m_values[d] = target[input->GetDim() +
                                     std::distance(names.begin(), field)];
            }
            else
            {
                std::vector<LibUtilities::FieldDefinitionsSharedPtr>
                    definitions;
                std::vector<std::vector<NekDouble>> data;
                LibUtilities::FieldIO::CreateForFile(m_session, filename)
                    ->Import(filename, definitions, data);
                std::vector<NekDouble> coefficients(m_expansion->GetNcoeffs(),
                                                    0.0);
                bool found = false;
                for (size_t f = 0; f < definitions.size(); ++f)
                    if (std::find(definitions[f]->m_fields.begin(),
                                  definitions[f]->m_fields.end(),
                                  variable) != definitions[f]->m_fields.end())
                    {
                        m_expansion->ExtractDataToCoeffs(
                            definitions[f], data[f], variable, coefficients);
                        found = true;
                    }
                ASSERTL0(found, "Missing magnetic field " + variable + " in " +
                                    filename);
                LibUtilities::Field<NekDouble, FieldState::Coeff> coeff(
                    "Magnetic file coefficients",
                    MultiRegions::GetBlockAttributes<NekDouble,
                                                     FieldState::Coeff>(
                        m_expansion),
                    1, 1);
                LibUtilities::Field<NekDouble, FieldState::Phys> phys(
                    "Magnetic file values",
                    MultiRegions::GetBlockAttributes<NekDouble,
                                                     FieldState::Phys>(
                        m_expansion),
                    1, 1);
                coeff.template CopyVector<NektarSpaces::HostSpace>(
                    coefficients);
                auto bwd = MultiRegions::BwdTransOp<NekDouble>::Create(
                    m_expansion, {name});
                bwd->Apply(coeff, phys);
                const auto values = phys.template ToVector<NekDouble>();
                std::copy(values.begin(), values.end(), m_values[d].data());
            }
        }
    }

    void ReadMeanFile()
    {
        LibUtilities::PtsFieldSharedPtr input;
        m_ptsIO.Import(m_session->GetFunctionFilename(m_functionName, "Bx"),
                       input);
        ASSERTL0(input->GetDim() == 2 && input->GetNFields() == 3,
                 "MagneticMeanField requires point data [R,Z,Br,BZ,Bphi].");
        Points points;
        input->GetPts(points);
        const size_t count = input->GetNpoints();
        ASSERTL0(count > 0, "MagneticMeanField point file is empty.");
        bool rebuild = !m_meanInput || m_meanCoordinates[0].size() != count;
        if (!rebuild)
            for (unsigned d = 0; d < 2; ++d)
                rebuild =
                    rebuild || !std::equal(points[d].begin(), points[d].end(),
                                           m_meanCoordinates[d].begin());
        if (rebuild)
        {
            m_meanCoordinates = Points(2);
            for (unsigned d = 0; d < 2; ++d)
            {
                m_meanCoordinates[d] = Array<OneD, NekDouble>(count);
                std::copy_n(points[d].data(), count,
                            m_meanCoordinates[d].data());
            }
            Points torus(6), target(6);
            for (unsigned d = 0; d < 6; ++d)
                torus[d] = Array<OneD, NekDouble>(m_increments * count, 0.0);
            for (unsigned d = 0; d < 3; ++d)
            {
                target[d]     = m_coords[d];
                target[d + 3] = m_values[d];
            }
            for (unsigned angle = 0; angle < m_increments; ++angle)
            {
                const double phi = angle * 2 * std::acos(-1.0) / m_increments;
                for (size_t p = 0; p < count; ++p)
                {
                    torus[0][angle * count + p] = std::cos(phi) * points[0][p];
                    torus[1][angle * count + p] = points[1][p];
                    torus[2][angle * count + p] = std::sin(phi) * points[0][p];
                }
            }
            m_meanInput =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                    3, std::vector<std::string>{"Bx", "By", "Bz"}, torus);
            m_meanOutput =
                MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                    3, std::vector<std::string>{"Bx", "By", "Bz"}, target);
            m_meanInterp.CalcWeights(m_meanInput, m_meanOutput);
        }
        Points torus;
        m_meanInput->GetPts(torus);
        for (unsigned angle = 0; angle < m_increments; ++angle)
        {
            const double phi = angle * 2 * std::acos(-1.0) / m_increments;
            for (size_t p = 0; p < count; ++p)
            {
                const size_t q = angle * count + p;
                torus[3][q] =
                    std::cos(phi) * points[2][p] - std::sin(phi) * points[4][p];
                torus[4][q] = points[3][p];
                torus[5][q] =
                    std::sin(phi) * points[2][p] + std::cos(phi) * points[4][p];
            }
        }
        m_meanInterp.Interpolate(m_meanInput, m_meanOutput);
    }
};

} // namespace PENKNIFE::detail
#endif
