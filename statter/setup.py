from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="mad_statter",
    version="0.1.0A",
    description="Read length statistics, alignment statistics etc...",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Sudeep Sahadevan",
    author_email="sahadeva@embl.de",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Environment :: Console",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
    install_requires=["bokeh==3.2.2", "pysam==0.21.0", "click==8.1.7", "pandas==2.1.0"],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "read_length_parser=scripts.read_length_stats:read_length_parser",
            "read_length_plotter=scripts.read_length_stats:read_length_plotter",
            "gene_type_read_length_parser=scripts.bam_stats:gene_type_read_length_stats",
            "gene_type_read_length_plotter=scripts.bam_stats:gene_type_length_plotter"
        ]
    },
)
