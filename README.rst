Biostructmap
============

Biostructmap is a Python tool for mapping sequence-aligned data (such as
location of polymorphisms) onto a protein structure.

Additionally, biostructmap allows for the incorporation of residue
spatial-proximity into sliding-window calculations, and can be used to
incorporate protein structure information into genetic tests of
selection pressure.

A web-based interface is available `here <https://biostructmap.burnet.edu.au>`__,
although the Python package is more flexible and likely to be faster.

Table of Contents
=================

-  `Usage Examples <#usage-examples>`__
-  `Prerequisites <#prerequisites>`__
-  `Installing <#installing>`__
-  `Testing <#running-the-tests>`__
-  `Contributing <#contributing>`__
-  `Versioning <#versioning>`__
-  `Authors <#authors>`__
-  `License <#license>`__
-  `Citing <#citing>`__
-  `Acknowledgments <#acknowledgments>`__

Getting Started
---------------

Usage Examples
--------------


Calculate proportion of polymorphic residues within a radius
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple usage case may be identification of regions of the protein with
a high percentage of polymorphic residues. If we are perhaps interested
in antibody-antigen interaction, 15 Angstrom is a reasonable radius over
which to average over.

::

    import biostructmap

    # Initialise structure object
    structure = biostructmap.Structure('1zrl.pdb', 'test_pdb_name')

    # The location of known polymorphisms relative to the PDB sequence (we are not
    # providing a reference sequence for this example), for each chain.
    data = {('A',): [200, 276, 300, 480, 367, 349]}

    # Map polymorphism data using a radius of 15 Angstrom. Results are returned
    # in a new object.
    results = structure.map(data, method='snps', ref=None, radius=15)

    # Use the results object to write data to a local PDB file, with data saved
    # in the B-factor column
    results.write_data_to_pdb_b_factor(fileobj='test_pdb_data_write.pdb')

Calculation of average hydrophobicity for all surface exposed residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A slighly more complicated usage case may be the calculation of an
average amino acid propensity scale, such as the Kyte & Doolittle index
of hydrophobicity. Additionally, if we are solely interested in surface
exposed residues, we may wish to restrict analysis to only residues with
a relative solvent accessibility greater than 0.2.

::

    import biostructmap

    # Initialise structure object
    structure = biostructmap.Structure('1zrl.pdb', 'test_pdb_name')

    # For this method, the data parameter is a string which represents the amino
    # acid propensity scale we wish to use. Note the use of the optional rsa_range
    # parameter to restrict to surface exposed residues.
    results = chain.map(data='kd', method='aa_scale', ref=None, radius=15,
                        rsa_range=(0.2, 1.0))

    # Use the results object to write data to a local PDB file, with data saved
    # in the B-factor column
    results.write_data_to_pdb_b_factor(fileobj='test_pdb_data_write.pdb')

Calculation of Tajima's D using protein structural information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can also use the biostructmap package to calculate a modified
Tajima's D value which incorporates protein structural information ---
essentially using a 3D sliding window instead of the standard 2D sliding
window often applied over a protein sequence.

::

    import biostructmap

    # Initialise structure object
    structure = biostructmap.Structure('1zrl.pdb', 'test_pdb_name')

    # Read in multiple sequence alignment data
    msa_data = biostructmap.SequenceAlignment('seq_align.fsa')
    data = {('A',): msa_data}

    # Reference seq might be the first sequence in the multiple sequence alignment
    reference_seq = {'A': str(msa_data[0].seq)}

    results = structure.map(data=data, method='tajimasd', ref=reference_seq,
                        radius=15, map_to_dna=True)

    results.write_data_to_pdb_b_factor(fileobj='test_pdb_data_write.pdb')

Result can be easily viewed in PyMol using the ``spectrum`` command.

From the Pymol command line:

::

    load my_pdb_file_name_here

    as surface

    #Select all residues with a mapped data value. Can change the default 'no-value'
    #option when writing to pdb b factor using biostructmap if needed.
    select nonzeros, b < 0 | b > 0

    color white

    spectrum b, selection=nonzeros

    #Make a publication quality image. May need to center molecule and perhaps
    #adjust image size to your requirements.
    set ray_opaque_background, off
    ray 2400, 2400
    cmd.png('output_file_name.png', dpi=300)

Prerequisites
-------------

Installing the biostructmap package requires both an install of the main
package, as well as optional install of a few external binaries (NCBI BLAST+,
Exonerate and DSSP).

BLAST+:
^^^^^^^

To install the BLAST+ package, visit the `NCBI BLAST+
site <https://blast.ncbi.nlm.nih.gov/>`__ and follow the links to
download and install a local copy of the BLAST+ application.

BLAST+ is not required, but is recommended. If BLAST+ is not installed,
a fallback pairwise alignment is performed using BioPython.pairwise2, and
the user should indicate that BLAST+ is not installed by including:

::

    import biostructmap

    biostructmap.seqtools.LOCAL_BLAST = False


DSSP:
^^^^^

To install DSSP, visit the `DSSP
website <http://swift.cmbi.ru.nl/gv/dssp/>`__ and follow the
instructions for install. Alternatively, users of recent Ubuntu or
Debian distributions will find that DSSP is available as part of these
distributions. To check if DSSP is currently installed under Linux, try
running:

::

    dssp --version || mkdssp --version

At least one of these should return version 2.x.x

If DSSP is not installed, you can try installing ``dssp`` using your
local package manager. For example, on Ubuntu:

::

    sudo apt-get install dssp

If this fails you will have to install DSSP from the source code
provided `here <http://swift.cmbi.ru.nl/gv/dssp/>`__.

DSPP is not strictly required, but any analysis that involves calculation
of secondary structure or solvent accessibility will raise an exception
if DSSP is not installed.

Exonerate:
^^^^^^^^^^

To install Exonerate, visit the `Exonerate
website <http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate>`__
and follow the instructions to install Exonerate on your system.
Alternatively, Exonerate is available through the default Ubuntu
repositories:

::

    sudo apt-get install exonerate

Note that Exonerate is only required if performing calculation of
Tajima's D over a protein structure using a multiple sequence alignment
- it is used to align a genomic sequence to a protein coding region. If
this functionality is not required, then biostructmap can be installed
and run without Exonerate, although some of the tests will fail.

If Exonerate is not installed, a fallback pairwise alignment is performed
using BioPython.pairwise2, and the user should indicate that Exonerate is not
installed by including:

::

    import biostructmap

    biostructmap.seqtools.LOCAL_EXONERATE = False

Numpy:
^^^^^^^^^^^^^

Before install biostructmap it is recommended to install Numpy
using your Python package manager of choice (eg pip or conda). If you
are using the Anaconda distribution of Python, then Numpy should be installed
already. If not, or if you are using a virtual environment:

::

    conda install numpy

or

::

    pip install numpy

SciPy:
^^^^^^^^^^^^^^

While there is no hard dependency on SciPy, calculation of nearby residues
can be very memory intensive without SciPy present. If you are getting a MemoryError
exception with large PDB files, then consider installing SciPy in your python environment.


Installing
----------

To install the biostructmap package, it is first recommended that you
make sure all tests pass in your environment.

From the root package directory, run:

::

    python setup.py test

If these tests pass, you can then install the package (or just skip
straight to this step if you're feeling lucky):

::

    python setup.py install

Running the tests
-----------------

From the root package directory run:

::

    python setup.py test

or alternatively

::

    pytest

These tests should cover most of the biostructmap functionality, with
several tests reliant on additional packages such as NCBI BLAST+ or
DSSP, which should be installed alongside biostructmap.

biostructmap was developed for Python 3+, but also supports Python 2.7.
Please contact us if any compatibility issues are observed with older
versions of Python.

Contributing
------------

Please read `CONTRIBUTING.rst <CONTRIBUTING.rst>`__ for details on our
code of conduct, and the process for submitting pull requests to us.

Versioning
----------

We use `SemVer <http://semver.org/>`__ for versioning. For the versions
available, see the `tags on this
repository <https://github.com/andrewguy/biostructmap/tags>`__.

Authors
-------

-  **Andrew Guy** - *Main Author* - `Github
   Page <https://github.com/andrewguy>`__

See also the list of
`contributors <https://github.com/andrewguy/biostructmap/contributors>`__
who participated in this project.

License
-------

This project is licensed under the MIT License - see the
`LICENSE.txt <LICENSE.txt>`__ file for details

Citing
------

If you have used this tool please cite:

-  Guy, A. J., Irani, V., Richards, J. S. & Ramsland, P. A. BioStructMap: A
   Python tool for integration of protein structure and sequence-based features.
   *Bioinformatics* (2018). doi:10.1093/bioinformatics/bty474

-  Guy, A. J. *et al.* Proteome-wide mapping of immune features onto
   Plasmodium protein three-dimensional structures. *Sci. Rep.* **8**, 4355 (2018).

Acknowledgments
---------------

-  Paul Ramsland, Jack Richards and Vashti Irani for various suggestions
   and support.
