from spack.package import *
import os
import shutil


class NesoTokamak(CMakePackage):
    """NESO-Tokamak"""

    homepage = "https://github.com/ExCALIBUR-NEPTUNE/NESO-Tokamak"

    git = "https://github.com/ExCALIBUR-NEPTUNE/NESO-Tokamak.git"
    version("working", branch="main")
    version("main", branch="main")

    depends_on("mpi", type=("build", "link", "run"))

    # Non-SYCL dependencies
    depends_on("hdf5 +mpi +hl", type=("build", "link", "run"))
    depends_on("cmake@3.24:", type="build")
    depends_on("cmake@3.24:", type="build", when="@main")
    depends_on("cmake@3.24:", type="build", when="@working")
    depends_on("neso.neso@working", type=("build", "link", "run"))
    depends_on("reactions", type=("build", "link", "run"))
    depends_on("neso.neso-rng-toolkit@working", type=("build", "link", "run"))

    # Depend on a sycl implementation - with workarounds for intel packaging.
    depends_on("sycl", type=("build", "link", "run"))
    depends_on("intel-oneapi-dpl", when="^dpcpp", type="link")
    conflicts("%dpcpp", msg="Use oneapi compilers instead of dpcpp driver.")
    conflicts("^dpcpp", when="%gcc", msg="DPC++ can only be used with Intel oneAPI compilers.")

    def cmake_args(self):
        args = []
        if not "+build_tests" in self.spec:
            args.append("-DENABLE_NESO_PARTICLES_TESTS=off")
        if "+nvcxx" in self.spec:
            args.append("-DNESO_PARTICLES_DEVICE_TYPE=GPU")
            args.append("-DHIPSYCL_TARGETS=cuda-nvcxx")

        return args
