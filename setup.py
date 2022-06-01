from setuptools import setup
import re

# taken from https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
VERSIONFILE = "icat/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name="icat",
    version=verstr,
    description="Identify cell-types across treatments in single-cell RNA sequencing data",
    url="https://github.com/dakota-hawkins/icat",
    author="Dakota Y. Hawkins",
    author_email="dyh0110@bu.edu",
    install_requires=[
        "ncfs",
        "sslouvain",
        "scikit-learn",
        "numpy",
        "scipy",
        "apricot",
        "pandas",
        "scanpy",
    ],
    license="BSD",
    packages=["icat"],
    zip_safe=False,
)
