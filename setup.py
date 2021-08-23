from setuptools import setup, find_packages
from codecs import open
from os import path
import beadpull.version
import sys


here = path.abspath(path.dirname(__file__))

if not sys.version_info[0] == 3 and sys.version_info[1] >= 5:
    sys.exit("Sorry, Python 3.5 or higher only")

# Get the long description from the README file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="beadpull",
    version=beadpull.version.version,
    author="S.V. Matsievskiy",
    author_email="matsievskiysv@gmail.com",
    maintainer="S.V. Matsievskiy",
    maintainer_email="matsievskiysv@gmail.com",
    url="https://gitlab.com/matsievskiysv/beadpull",
    description="Bead pulling measurement software for RF cavity measurements",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GPLv3+",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: X11 Applications :: Qt",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public \
License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="accelerators particles measurements",
    packages=find_packages(
        exclude=["*.tests", "*.tests.*", "tests.*", "tests"]
    ),
    include_package_data=True,
    install_requires=["numpy", "PyQt5", "PyQtWebEngine", "matplotlib", "pysmithplot-fork", "PyYAML", "PyVISA", "pyvisa-py"],
    extras_require={"dev": ["check-manifest"]},
    entry_points={"console_scripts": ["beadpull = beadpull.__main__:main"]},
)
