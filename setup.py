"""Installation script."""

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

# ## Provides python function "numpy.get_include()" that returns
# ## the path to the core numpy header file,
# ## which is required for compilation.
import numpy

# ## Required for cython translation.
# ## Otherwise search for file "eigency/core.pxd" fails.
import sys
sys.path.append("cy/eigency")

# ## Provides python function "eigency.get_include()" which returns
# ## the path to various headers, such as "eigency_cpp.h",
# ## which are required for compilation.
import eigency

extensions = [
    Extension(
        "hntf2d.map",
        sources=["cy/hntf2d/map.pyx", "cpp/src/map.cpp"],
        language="c++",
        extra_compile_args=["-O2", "-std=c++11"],
        include_dirs=[
            "cpp/include",
            "cy/hntf2d",
            numpy.get_include()
        ] + eigency.get_includes(),
        libraries=[
            "tinyxml2"
        ]
    ),
    Extension("hntf2d.utils", ["cy/hntf2d/utils.py"],
              language="python")
]

setup(
    name="hntf2d",
    version="0.1.dev1",
    ext_modules=cythonize(extensions,
                          gdb_debug=False,
                          compiler_directives={"language_level": 3}),
    packages={"cy/hntf2d": "hntf2d"},
    cmdclass={build_ext: build_ext},
    install_requires=['eigency'],
    dependency_links=['cy/eigency']
)
