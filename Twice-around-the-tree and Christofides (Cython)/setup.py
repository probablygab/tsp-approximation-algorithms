# cython: boundscheck=False, initializedcheck=False, nonecheck=False, cdivision=True
from setuptools import setup
from Cython.Build import cythonize
from Cython.Compiler import Options

Options.gcc_branch_hints = True
Options.convert_range = True

setup(
    ext_modules = cythonize("tp2.py", annotate=True)
)