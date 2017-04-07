# StructMap

Structmap is a Python tool for mapping sequence-aligned data (such as location of polymorphisms) onto a protein structure.

Additionally, Structmap allows for the incorporation of residue spatial-proximity into sliding-window calculations, and can be used to incorporate protein structure information into genetic tests of selection pressure.

## Getting Started


## Usage Examples

A few usage examples:

```
import structmap

structure = structmap.Structure('1zrl.pdb', 'test_pdb_name')

data = [200, 276, 300, 480, 367, 349]

chain = structure[0]['A']

results = chain.map(data, method='snps', ref=None, radius=15)

results.write_data_to_pdb_b_factor(filename='test_pdb_data_write.pdb')
```

## Prerequisites

Installing the Structmap package requires both an install of the main package, as well as install of a few external binaries (NCBI BLAST+, Exonerate and DSSP).

#### BLAST+:

To install the BLAST+ package, visit the [NCBI BLAST+ site](https://blast.ncbi.nlm.nih.gov/) and follow the links to download and install a local copy of the BLAST+ application.

#### DSSP:

To install DSSP, visit the [DSSP website](http://swift.cmbi.ru.nl/gv/dssp/) and follow the instructions for install. Alternatively, users of recent Ubuntu or Debian distributions will find that DSSP is available as part of these distributions. To check if DSSP is currently installed under Linux, try running:

```
dssp --version || mkdssp --version
```

At least one of these should return version 2.x.x

If DSSP is not installed, you can try installing `dssp` using your local package manager. For example, on Ubuntu:

```
sudo apt-get install dssp
```

If this fails you will have to install DSSP from the source code provided [here](http://swift.cmbi.ru.nl/gv/dssp/).

#### Exonerate:

To install Exonerate, visit the [Exonerate website](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) and follow the instructions to install Exonerate on your system. Alternatively, Exonerate is available through the default Ubuntu repositories:

```
sudo apt-get install exonerate
```

Note that Exonerate is only required if performing calculation of Tajima's D over a protein structure using a multiple sequence alignment - it is used to align a genomic sequence to a protein coding region. If this functionality is not required, then Structmap can be installed and run without Exonerate, although some of the tests will fail.

#### Numpy, Scipy:

Before install Structmap it is recommended to install Numpy and Scipy using your Python package manager of choice (eg pip or conda). If you are using the Anaconda distribution of Python, then both NumPy and SciPy should be installed already. If not, or if you are using a virtual environment:

```
conda install numpy scipy
```

or

```
pip install numpy scipy
```

## Installing

To install the Structmap package, it is first recommended that you make sure all tests pass in your environment.

From the root package directory, run:

```
python setup.py test
```

If these tests pass, you can then install the package (or just skip straight to this step if you're feeling lucky):

```
python setup.py install
```

## Running the tests

From the root package directory run:

```
python setup.py test
```

or alternatively

```
pytest
```

These tests should cover most of the structmap functionality, with several tests reliant on additional packages such as NCBI BLAST+ or DSSP, which should be installed alongside structmap.

Structmap was developed for Python 3+, but also supports Python 2.7. Please contact us if any compatibility issues are observed with older versions of Python.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/andrewguy/structmap/tags).

## Authors

* **Andrew Guy** - *Main Author* - [Github Page](https://github.com/andrewguy)

See also the list of [contributors](https://github.com/andrewguy/structmap/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Paul Ramsland, Jack Richards and Vashti Irani for various suggestions and support.
