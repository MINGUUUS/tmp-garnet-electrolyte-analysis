"""
Setup script for garnet-electrolyte-analysis package.
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read version from __init__.py
version = {}
with open("garnet_analysis/__init__.py") as fp:
    exec(fp.read(), version)

setup(
    name="garnet-electrolyte-analysis",
    version=version["__version__"],
    author=version["__author__"],
    author_email=version["__email__"],
    description="A comprehensive Python toolkit for computational analysis of garnet-structured solid electrolytes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/garnet-electrolyte-analysis",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "matplotlib>=3.5.0",
        "tqdm>=4.62.0",
    ],
    extras_require={
        "full": [
            "pymatgen>=2022.0.0",
            "ase>=3.22.0",
        ],
        "aimd": [
            "aimd",  # If available through PyPI
        ],
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=22.0",
            "flake8>=4.0",
            "mypy>=0.910",
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
        "notebooks": [
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",
            "ipywidgets>=7.6.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "garnet-analyze=garnet_analysis.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "garnet_analysis": [
            "data/*.csv",
            "data/*.cif",
            "examples/*.ipynb",
        ]
    },
    keywords="materials science, solid electrolytes, garnet, diffusivity, molecular dynamics, battery materials",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/garnet-electrolyte-analysis/issues",
        "Documentation": "https://garnet-electrolyte-analysis.readthedocs.io/",
        "Source": "https://github.com/yourusername/garnet-electrolyte-analysis",
    },
)