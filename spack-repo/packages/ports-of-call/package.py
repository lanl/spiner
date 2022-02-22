# Spackage for Ports of Call

from spack import *

class PortsOfCall(CMakePackage, CudaPackage):
    """ports-of-call"""

    homepage    = "https://github.com/lanl/ports-of-call"
    url         = "https://github.com/lanl/ports-of-call/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/ports-of-call.git"

    version("main", branch="main")
    variant("doc", default=False, description="Sphinx Documentation Support")
    variant("portability_strategy", description="Portability strategy backend",
            values=("Kokkos","Cuda","None"),multi=False,default="kokkos",
            when="^kokkos")
    variant("portability_strategy", description="Portability strategy backend",
            values=("Kokkos","Cuda","None"),multi=False,default="none")

    depends_on("cmake@3.12:")

    depends_on("py-sphinx", when="+doc")
    depends_on("py-sphinx-rtd-theme@0.4.3", when="+doc")
    depends_on("py-sphinx-multiversion", when="+doc")

    def cmake_args(self):
        args = [
            self.define_from_variant("PORTABILITY_STRATEGY", "portability_strategy")
        ]
        return args
