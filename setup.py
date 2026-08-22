#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import sys

# --- Extract version from version.py without loading the package ---
# This is the most reliable way if the file exists
try:
    version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_mut_integrate', 'version.py')
    with open(version_file) as f:
        exec(f.read())  # defines __version__
except (FileNotFoundError, NameError):
    # Fallback version if version.py is missing
    __version__ = "0.1.0"

# --- Read the long description from README.md ---
try:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "Integrate wild-type and mutant EnzyWizard reports into graph-based structures."

# --- Setup configuration ---
setup(
    name="enzywizard-mut-integrate",
    version=__version__,
    author="bioinfbrad",
    description=(
        "Integrate wild-type and mutant EnzyWizard JSON reports to build structured paired "
        "graph datasets, enabling mutation effect analysis and graph-based machine learning."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-mut-integrate",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    # Core runtime dependencies based on the tool's functionality
    install_requires=[
        "biopython>=1.86",          # Sequence handling, residue mapping
        "numpy>=1.23.5,<2",         # Numerical operations
        "packaging",                # Version handling
        # 'rdkit', 'openmm', 'prody', 'fair-esm', 'vina' are NOT required here,
        # because they are run-time dependencies of the individual analysis tools.
        # This tool only integrates their JSON outputs.
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-mut-integrate = enzywizard_mut_integrate.cli:main",
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
