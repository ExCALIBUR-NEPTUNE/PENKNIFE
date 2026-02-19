#ifndef MAGNETIC_FIELD_HPP
#define MAGNETIC_FIELD_HPP

#include <FieldUtils/Interpolator.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Core/SessionFunction.h>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;

namespace PENKNIFE
{
class PlasmaSystem;

class MagneticField
{
    friend PlasmaSystem;

public:
    MagneticField(const LU::SessionReaderSharedPtr &session,
                  const std::weak_ptr<PlasmaSystem> &eq_sys,
                  Array<OneD, MR::DisContFieldSharedPtr> &B, const int dim);

    void Read(NekDouble time = 0);
    void Solve(Array<OneD, Array<OneD, NekDouble>> &J = NullNekDoubleArrayOfArray);
    void Update(NekDouble time = 0);

private:
    static const inline std::vector<std::string> Bstring = {"Bx", "By", "Bz"};
    int ndim;
    int npoints;

    enum class field_type : int
    {
        mean_file,
        mean_exp,
        func,
        solve
    } type;

    SU::SessionFunctionSharedPtr m_function;

    // For mean_exp
    LibUtilities::EquationSharedPtr Bxfunc;
    LibUtilities::EquationSharedPtr Byfunc;
    LibUtilities::EquationSharedPtr Bzfunc;
    Array<OneD, NekDouble> r;
    Array<OneD, NekDouble> phi;
    Array<OneD, NekDouble> x1;
    Array<OneD, NekDouble> Br;
    Array<OneD, NekDouble> Bphi;

    //for mean file
    FieldUtils::Interpolator<std::vector<MR::ExpListSharedPtr>> interp;
    LU::PtsIO ptsIO;
    std::string filename;
    Array<OneD, Array<OneD, NekDouble>> B2D;
    Array<OneD, Array<OneD, NekDouble>> B3D;
    LU::PtsFieldSharedPtr inPts;
    LU::PtsFieldSharedPtr outPts;
    unsigned int increments = 64;

    bool transient_field;
    /// Current density
    Array<OneD, MR::DisContFieldSharedPtr> J;

    /// Magnetic field vector
    Array<OneD, MR::DisContFieldSharedPtr> B;
    Array<OneD, Array<OneD, NekDouble>> B_in;

    /// Normalised magnetic field vector
    Array<OneD, Array<OneD, NekDouble>> b_unit;

    /// Squared Magnitude of the magnetic field
    Array<OneD, NekDouble> mag_B;

    std::weak_ptr<PlasmaSystem> m_sys;
};
} // namespace PENKNIFE

#endif