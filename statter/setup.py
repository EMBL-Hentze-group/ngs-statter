from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name = 'mad_statter',
    version = '0.1.0A',
    description = 'Read length statistics, alignment statistics etc...',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author = 'Sudeep Sahadevan',
    author_email = 'sahadeva@embl.de',
    license = 'MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Environment :: Console',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    #install_requires = ['pytorch==1.13.0', 'pytorch-cuda==11.7', 'ray-all==2.1.0', ' pyyaml==6.0', 'sentencepiece==0.1.95', 'transformers==4.24.0', 'scikit-learn==1.1.3'],
    packages=find_packages(),
    entry_points = {
        'console_scripts': [
            'read_length_paser=scripts.parse_read_length:read_length_parser'
            ]
    }
)
