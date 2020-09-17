"""
Parameters for module setup
"""

from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name = "Fibroblast_Gene_Regulatory_Network",
    version = 2020.09,
    description = "Inference, validation, and analysis of transcriptional regulatory network derived from myofibroblast RNA-seq data",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = "https://github.com/jdroger/Fibroblast_Gene_Regulatory_Network",
    author = "Jesse Rogers",

    classifiers = [
        "Development Status :: 3 - Alpha"
        "Programming Language :: Python :: 3.8"
    ],

    package_dir={'': 'src'},
    python_requires='>=3.7, <4',
    install_requires=['numpy', 'pandas', 'dask', 'distributed', 'arboreto']

)