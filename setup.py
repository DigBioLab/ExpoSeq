from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
#module = Extension ('loop_morosita_horn', sources=[r'C:\Users\nils\PycharmProjects\ExpoSeq\python_scripts\plots\create_heatmap.pyx'])

setup(
    ext_modules = cythonize('create_heatmap.pyx'),
    include_dirs=[numpy.get_include()],
)

#[build - system]
#requires = ["setuptools", "wheel", "Cython"]