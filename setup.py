from setuptools import setup, find_packages
from codecs import open
from pathlib import Path
import bi2d.version


here = Path(__file__).parent.absolute()

# Get the long description from the README file
with open(here / "README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="BeamImpedance2DX",
    version=bi2d.version.version,
    author="S.V. Matsievskiy",
    author_email="matsievskiysv@gmail.com",
    maintainer="S.V. Matsievskiy",
    maintainer_email="matsievskiysv@gmail.com",
    url="https://gitlab.com/matsievskiysv/bi2d",
    description="Calculate electromagnetic field and beam impedance in 2D device section",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="AGPLv3+",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="accelerators particles impedance electromagnetic",
    packages=find_packages(
        exclude=["*.tests", "*.tests.*", "tests.*", "tests"]
    ),
    include_package_data=True,
    install_requires=["fenics-dolfinx>=0.3.1", "python_version>=3.6"],
    extras_require={"dev": ["check-manifest"]},
)
