from __future__ import absolute_import, division, print_function

import io
import os
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

# Version is defined in structmap/__init__.py
__version__ = "Undefined"
for line in open('structmap/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

LONG_DESCRIPTION = read('README.txt', 'CHANGES.txt')

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
    version=__version__,
    url='',
    author='Andrew Guy',
    tests_require=['pytest'],
    setup_requires=['numpy'],
    install_requires=['Biopython>=1.66',
                      'DendroPy>=4.0.3',
                      'numpy',
                      'scipy'
                     ],
    cmdclass={'test': PyTest},
    author_email='andrewguy@burnet.edu.au',
    description='A simple package for mapping data onto protein PDB structures',
    long_description=LONG_DESCRIPTION,
    packages=['structmap'],
    include_package_data=True,
    platforms='any',
    test_suite='tests.test_structmap',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
