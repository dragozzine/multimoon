from distutils.core import setup
from Cython.Build import cythonize
import sysconfig
import os

os.environ['CFLAGS'] = '-Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -std=c++11'

setup(
        name = "spinny",
        ext_modules = cythonize(('spinny.pyx')),
        extra_compile_args=["-std=c++11"],
        language='c++11'
        )