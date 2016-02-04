from __future__ import absolute_import, division, print_function

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

import structmap

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='structmap',
    version=structmap.__version__,
    url='',
    license='Apache Software License',
    author='Andrew Guy',
    tests_require=['pytest'],
    install_requires=['Biopython>=1.66',
                    'DendroPy>=4.0',
                    'NumPy',
                    'SciPy'
                    ],
    cmdclass={'test': PyTest},
    author_email='andrewguy@burnet.edu.au',
    description='A simple package for mapping data onto protein PDB structures',
    long_description=long_description,
    packages=['structmap'],
    include_package_data=True,
    platforms='any',
    test_suite='tests.test_structmap',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
