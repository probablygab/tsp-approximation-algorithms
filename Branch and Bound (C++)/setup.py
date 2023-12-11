import pybind11
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Compile the C++ code into a Python module
setup(ext_modules = [
        Extension(
            'bnb',                                  # Module name in Python
            ['BranchAndBound.cpp'],                 # List of C++ source files
            include_dirs=[pybind11.get_include()],  # Include pybind11 headers
            language='c++'
        )
    ]
)