import os
import pathlib
import platform
import re
import subprocess
import sys
import sysconfig

from packaging.version import Version
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

here = pathlib.Path(__file__).parent.resolve()
# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

EXTENSION_CFG = {"define_macros": []}
if sys.version_info >= (3, 12):
    EXTENSION_CFG["py_limited_api"] = True
    EXTENSION_CFG["define_macros"].append(
        ("Py_LIMITED_API", 0x030C0000 if sys.version_info < (3, 15) else 0x030F0000)  # PY_VERSION_HEX
    )


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[], **EXTENSION_CFG)
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError as err:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from err

        if platform.system() == "Windows":
            cmake_version = Version(re.search(r"version\s*([\d.]+)", out.decode()).group(1))
            if cmake_version < Version("3.1.0"):
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}",
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DPYBIND11_STABLE_ABI=ON",
            "-DPYBIND11_ABI3T=ON",
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE:STRING=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=self.build_temp)
        subprocess.check_call(["cmake", "--install", "."], cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


SETUPTOOLS_OPTIONS = {}
if sys.version_info >= (3, 12):
    SETUPTOOLS_OPTIONS["bdist_wheel"] = {"py_limited_api": "cp312" if sys.version_info < (3, 15) else "cp315"}

setup(
    name="lephare",
    url="https://lephare.readthedocs.io/en/latest/",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # add extension module
    ext_modules=[CMakeExtension("lephare._lephare")],
    include_package_data=True,
    # add custom build_ext command
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    python_requires=">=3.9",
    extras_require={"test": "pytest"},
    options=SETUPTOOLS_OPTIONS,
)
