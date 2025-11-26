# coding-utf-8

import os
from setuptools import setup, find_packages


if os.path.exists('README.md'):
    long_description = open('README.md').read()
else:
    long_description = '''ARKP-DB: An Automated Database of Astrochemical Reactions and Kinetics Parameters extracted from scientific literature.'''

setup(
    name='astrochemreactionextractor',
    version='1.0',
    license='MIT',
    packages=find_packages(),
    description='ARKP-DB: An Automated Database of Astrochemical Reactions and Kinetics Parameters extracted from scientific literature',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='text-mining mining Astrochemical Reactions Kinetics Parameters nlp science scientific',
    author='',
    author_email='',
    url='https://github.com/***',
    zip_safe=False,
    install_requires=[
        'PyMuPDF==1.26.1',
        'jellyfish==1.2.1',
        'matplotlib==3.10.7',
        'numpy==2.3.5',
        'pandas==2.3.3',
        'rdkit==2025.3.6',
        'requests==2.32.4',
        'sympy==1.14.0',
        'mpmath==1.3.0',
        'tqdm==4.67.1'
    ],
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Astrochemical',
        'Topic :: Text Processing',
    ],
)