from setuptools import setup
import os
import sys
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "vendor", "pybind11"))
from pybind11.setup_helpers import Pybind11Extension
del sys.path[-1]


__version__ = "0.0.1"
prefix = os.getenv('PREFIX')
source = os.getenv('LEPINJSOURCEPATH')

ext_modules = [
    Pybind11Extension("Geometry",
        [source+"projects/geometry/private/pybindings/geometry.cxx"],
        define_macros = [('VERSION_INFO', __version__)],
        library_dirs=[prefix+'/lib/'],
        libraries=['LeptonInjector'],
        ),
    Pybind11Extension("EarthModel",
        [source+"projects/detector/private/pybindings/earthmodel.cxx"],
        define_macros = [('VERSION_INFO', __version__)],
        library_dirs=[prefix+'/lib/'],
        libraries=['LeptonInjector'],
        ),
]

setup(
    name="LeptonInjector_PythonBindings",
    version=__version__,
    author="Nicholas Kamp",
    author_email="nwkamp@mit.edu",
    url="https://github.com/austinschneider/LeptonInjectorDevPrivate",
    description="pybind11 python bindings for LeptonInjector",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    zip_safe=False,
    python_requires=">=3.7",
)
