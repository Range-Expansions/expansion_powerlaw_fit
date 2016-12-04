from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import cython_gsl
import numpy as np

extensions = [
    Extension("expansion_powerlaw_fit",
              sources=["expansion_powerlaw_fit.pyx"],
              language="c", libraries = cython_gsl.get_libraries(),
              library_dirs = [cython_gsl.get_library_dir()],
              include_dirs = [cython_gsl.get_cython_include_dir(), np.get_include()])
]

setup(
    name='expansion_powerlaw_fit',
    version='1.0',
    url='',
    license='',
    author='Bryan Weinstein',
    author_email='bweinstein@seas.harvard.edu',
    description='Techniques to extract powerlaws written in Cython.',
    include_dirs = [cython_gsl.get_include(), np.get_include()],
    ext_modules = cythonize(extensions, annotate=True)
)