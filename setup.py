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

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="icat-sc",
    version=verstr,
    description="Identify cell states across treatments in single-cell RNA sequencing experiments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BradhamLab/icat",
    author="Dakota Y. Hawkins",
    author_email="dyh0110@bu.edu",
    install_requires=[
        "apricot-select",
        "igraph",
        "ncfs",
        "numpy",
        "pandas",
        "scanpy",
        "scikit-learn",
        "scipy",
        "louvain",
        "sslouvain",
    ],
    license="BSD",
    packages=["icat"],
    zip_safe=False,
)
