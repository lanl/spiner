# Spackage for Spiner Dependencies

from spack import *

class SpinerDeps(BundlePackage):
    """Spiner Dependencies"""

    homepage    = "https://github.com/lanl/spiner"
    url         = "https://github.com/lanl/spiner/archive/refs/heads/main.zip"
    git         = "git@github.com:lanl/spiner.git"

    version("main", branch="main")

    variant('doc', default=False, description='Enable Sphinx Documentation Support')

    depends_on("cmake@3.12:")
    depends_on("catch2@2.13.4:2.13.6")

    depends_on('py-sphinx', when='+doc')
    depends_on('py-sphinx-rtd-theme@0.4.3', when='+doc')
    depends_on('py-sphinx-multiversion', when='+doc')

