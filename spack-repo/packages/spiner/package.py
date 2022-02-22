# Spackage for Spiner

from spack import *

class Spiner(CMakePackage, CudaPackage):
    """Spiner"""

    homepage    = "https://github.com/lanl/spiner"
    url         = "https://github.com/lanl/spiner/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/spiner.git"

    version("main", branch="main")

    variant("kokkos", default=True, description="Enable kokkos",
            when="^kokkos")
    variant("openmp", default=True, description="Enable openmp kokkos backend",
            when="^openmp")
    variant("cuda", default=False, description="Enable cuda backend",
            when="^cuda")
    variant("kokkos", default=False, description="Enable kokkos",)
    variant("openmp", default=False, description="Enable openmp kokkos backend")
    variant("cuda", default=False, description="Enable cuda backend")

    variant("python", default=False, description="Python, Numpy & Matplotlib Support")
    variant("doc", default=False, description="Sphinx Documentation Support")
    variant("format", default=False, description="Clang-Format Support")

    depends_on("cmake@3.12:")
    depends_on("catch2@2.13.4:2.13.6")

    # Currently the raw cuda backend of ports-of-call is not supported.
    depends_on("ports-of-call portability_strategy=Kokkos", when="+kokkos")
    depends_on("ports-of-call portability_strategy=None", when="~kokkos")
    for _flag in list(CudaPackage.cuda_arch_values):
        depends_on("kokkos@3.2.00 cuda_arch=" + _flag, when="+cuda+kokkos cuda_arch=" + _flag)
    for _flag in ("~cuda", "+cuda", "~openmp", "+openmp"):
        depends_on("kokkos@3.2.00" + _flag, when="+kokkos" + _flag)
    depends_on("kokkos@3.2.00~shared+wrapper+cuda_lambda+cuda_relocatable_device_code", when="+cuda+kokkos")

    depends_on("python", when="+python")
    depends_on("py-numpy", when="+python")
    depends_on("py-matplotlib", when="+python")

    depends_on("py-sphinx", when="+doc")
    depends_on("py-sphinx-rtd-theme@0.4.3", when="+doc")
    depends_on("py-sphinx-multiversion", when="+doc")

    depends_on("llvm@12.0.0+clang", when="+format")

    conflicts("+cuda", when="~kokkos")
    conflicts("+openmp", when="~kokkos")
    conflicts("cuda_arch=none", when="+cuda", msg="CUDA architecture is required")

    def cmake_args(self):
        args = [
            "-DBUILD_TESTING={0}".format("ON" if self.run_tests else "OFF"),
        ]
        return args
