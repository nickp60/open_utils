"""
Setup for open_utils
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
import re
from codecs import open
from os import path
import sys

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))


if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: clermontpcr requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

## parse requirements file
install_reqs = parse_requirements("requirements.txt",
                                  session=False)
requirements = [str(ir.req) for ir in install_reqs]

setup(
    name='open_utils',
    version="0.0.3",

    description="you didn't think you needed it til you did",
    # long_description=long_description,
    long_description="""
    check out the GitHub
    repo for the real README.md file
    """,

    url='https://github.com/nickp60/EzClermont',

    # Author details
    author='Nick Waters',
    author_email='nickp60@gmail.com',
    license='MIT',
    # handle requirments
    install_requires=requirements,
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics evolution genomics development',
    packages=find_packages(),
    # add stuff to the MANIFEST
    include_package_data=True,
    entry_points={
       'console_scripts': [
           'extractRegion=extractRegion.extractRegion:main',
           'snagnblast=virulenceParser.snagnblast.snagnblast:main'
       ],
    },
)
