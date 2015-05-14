T\ :sub:`R`\ ENTo
=================
*Reduced Thickness Event-by-event Nuclear Topology*

.. image:: http://img.shields.io/travis/Duke-QCD/trento.svg?style=flat-square
  :target: https://travis-ci.org/Duke-QCD/trento

T\ :sub:`R`\ ENTo is a simple, fast model for the initial conditions of high-energy nuclear collisions (pp, pA, AA).
`arXiv:1412.4708 [nucl-th] <http://inspirehep.net/record/1334386>`_ formally presents the model and preliminary results.

Installation
------------
Prerequisites:

- `CMake <http://www.cmake.org>`_ 2.8+
- A `C++11 compiler <http://en.cppreference.com/w/cpp/compiler_support>`_ (preferably GCC 4.8+ or Clang 3.3+)
- The `Boost <http://www.boost.org>`_ C++ libraries, including runtime components  ``filesystem`` and ``program_options``
- (optional) The `HDF5 <http://www.hdfgroup.org/HDF5>`_ C++ library

All these dependencies are readily available on any operating system.
Some example installation commands for a few Linux distributions:

Ubuntu::

   apt-get install cmake g++ libboost-dev libboost-{filesystem,program-options}-dev libhdf5-dev

Fedora::

   yum install cmake gcc-c++ boost boost-{filesystem,program_options} hdf5 hd5-devel

Arch::

   pacman -S cmake gcc boost boost-libs hdf5-cpp-fortran

After installing the dependencies, download the `latest release <https://github.com/Duke-QCD/trento/releases/latest>`_ or clone the repository, then compile and install T\ :sub:`R`\ ENTo through the standard CMake sequence::

   mkdir build && cd build
   cmake ..
   make install

This will install the compiled binary to ``~/.local/bin/trento``.
If you do not want this to happen, run ``make`` instead of ``make install`` and the binary will be left at ``build/src/trento``.
The remainder of this document assumes ``trento`` is in your ``PATH``.

The code is `continuously tested <https://travis-ci.org/Duke-QCD/trento>`_ on Ubuntu with GCC 4.9.
It should run just as well on any Linux distribution or OS X, and probably on Windows.
For the compiler, Clang works as well as GCC.
Other compilers should work but may require modifying the compiler flags.

Usage
-----
T\ :sub:`R`\ ENTo has a standard command-line interface.
The basic syntax is ::

   trento [options] projectile projectile [number-events = 1]

where the only required arguments are the two projectile names.
For example, ``trento Pb Pb 10`` would run ten lead-lead events.

The remaining optional arguments may be given in any order, before or after the projectiles.
Run ``trento --help`` for a brief summary of the options and see below for more detailed descriptions and some `Examples`_.

Specifying projectiles
~~~~~~~~~~~~~~~~~~~~~~
The ``projectile`` arguments take species abbreviations, e.g. ``p``, ``Pb``, etc.
The known species are

======  =======  ============  ========
Symbol  Name     No. nucleons  Deformed
======  =======  ============  ========
p       proton   1             ---
Cu      copper   62            no
Cu2     copper   62            yes
Au      gold     197           no
Au2     gold     197           yes
Pb      lead     208           no
U       uranium  238           yes
======  =======  ============  ========

All species except the proton sample nucleons from a `Woods-Saxon <https://en.wikipedia.org/wiki/Woods%E2%80%93Saxon_potential>`_ distribution, either spherically symmetric or deformed as indicated.
The naming convention and Woods-Saxon parameters follow the `PHOBOS Glauber <http://inspirehep.net/record/1310629>`_ model.

General options
~~~~~~~~~~~~~~~
These are general options that don't fit in any other category.

-h, --help
   Show the help message and exit.

--version
   Print version number and exit.

--bibtex
   Print bibtex entry and exit.

-c, --config-file FILE
   Path to configuration file (see `Configuration files`_ below).
   May be given multiple times.


Output options
~~~~~~~~~~~~~~
The default output mode is to print event-by-event properties to stdout, in the following order::

   event_number impact_param npart mult e2 e3 e4 e5

with one line for each event.
``mult`` is the total initial entropy and the ``en`` are eccentricity harmonics ɛ\ :sub:`n`.
This format is designed for easy parsing, redirection to files, etc.

By default, the actual initial entropy profiles (grids) are not output.
There are two available output formats: text and HDF5 (if compiled).

In text mode, each event is written to a separate text file.
Each file has a commented header containing the event properties, like this::

   # event 0
   # b     = 2.964077155
   # npart = 380
   # mult  = 168.603282
   # e2    = 0.01953253866
   # e3    = 0.08961920965
   # e4    = 0.1101683349
   # e5    = 0.1727159106

The profile follows the header as a standard block-style grid.

HDF5 is a high-performance, cross-platform binary format for large numerical datasets.
Libraries are available in `most languages <https://en.wikipedia.org/wiki/Hierarchical_Data_Format#Interfaces>`_.
HDF5 is significantly faster than text output:
writing an event to a text file usually takes much longer than computing the actual event;
writing to HDF5 incurs only a small overhead.
Therefore, HDF5 is the recommended output format.

In HDF5 mode, all events are written to a single file with each event in a separate HDF5 dataset.
Event properties are written to each dataset as HDF5 attributes with names ``b``, ``npart``, ``mult``, ``e2``, etc.

-q, --quiet
   Disable printing event properties to stdout.
   Since both text and HDF5 output contain the event properties, it's often desirable to specify this option along with the output option.

-o, --output PATH
   Path to output events.
   If the path has an HDF5-like extension (``.hdf5``, ``.hdf``, ``.hd5``, ``.h5``), then all events will be written to that HDF5 file.
   Otherwise, the path is interpreted as a directory and events will be written to numbered text files in the directory.

   For text output, the directory will be created if it does not exist.
   If it does already exist, it must be empty (this is to avoid accidentally overwriting files or spewing thousands of files into an already-used location).

   For HDF5 output, the file must not already exist.
   Each event will be written as a numbered dataset in the file, and the standard event properties will be written as dataset attributes.

   Example:

   - ``--output events`` will write to text files ``events/0.dat``, ``events/1.dat``, ...
   - ``--output events.hdf`` will write to HDF5 file ``events.hdf`` with dataset names ``event_0``, ``event_1``, ...

Physical options
~~~~~~~~~~~~~~~~
These options control the physical behavior of the model.
They all have reasonable defaults, however **the defaults are not in any way a best-fit to experimental data**.
They are simply round numbers.
It is entirely expected that the ideal parameters will change depending on the beam energy.
In particular, **the cross section must be explicitly set for each beam energy**.

-p, --reduced-thickness FLOAT
   Reduced thickness parameter *p*.
   The reduced thickness is defined as the `generalized mean <https://en.wikipedia.org/wiki/Generalized_mean>`_ of participant nuclear thickness

   .. image:: http://latex2png.com/output//latex_11011000a8160e4838e75a0c11f293b2.png

   The default is *p* = 0, which corresponds to the geometric mean.

-k, --fluctuation FLOAT
   `Gamma distribution <https://en.wikipedia.org/wiki/Gamma_distribution>`_ shape parameter *k* for nucleon fluctuations.
   Fluctuations are sampled from a gamma distribution with the scale parameter fixed so that the mean is one:

   .. image:: http://latex2png.com/output//latex_17f24b3c97fb2b649d3dc4de4cd7e026.png

   The default is *k* = 1, which corresponds to an exponential distribution.
   For small *k*, the distribution has a long tail, leading to large fluctuations.
   For large *k*, the distribution becomes a narrow Gaussian, and eventually a delta function for very large values.

-w, --nucleon-width FLOAT
   Gaussian nucleon width in fm:

   .. image:: http://latex2png.com/output//latex_0c9ba0458eb84402a2a0fe505dc7164d.png

   The default is 0.5 fm.
   A reasonable range is roughly 0.4–0.8 fm.

-n, --normalization FLOAT
   Overall normalization factor.
   The default is 1.

-x, --cross-section FLOAT
   Inelastic nucleon-nucleon cross section σ\ :sub:`NN` in fm\ :sup:`2`.
   The default is 6.4 fm\ :sup:`2`, which is the approximate experimental value at LHC energy, √s = 2.76 TeV.

--b-min FLOAT
   Minimum impact parameter.
   The default is zero.

--b-max FLOAT
   Maximum impact parameter.
   The default is to run minimum-bias collisions for the given collision system.

   To run at fixed impact parameter, give the same value for both the min and the max.

--random-seed POSITIVE_INT
   Primarily for testing and debugging.

Grid options
~~~~~~~~~~~~
The thickness functions are discretized onto a square *N* × *N* grid centered at (0, 0).
The grid can have a dramatic effect on code speed and precision, so should be set carefully.
Computation time is roughly proportional to the number of grid cells (i.e. *N*\ :sup:`2`).

--grid-max FLOAT
   *x* and *y* maximum of the grid in fm, i.e. the grid extends from -max to +max.
   The default is 10 fm, large enough to accommodate all collision systems.
   However, this should be set as small as possible, since an unnecessarily large grid slows down the code.
   For anything but uranium-uranium, 9 fm is sufficient.
   For pp and pA, 3 fm is usually a good choice.

--grid-step FLOAT
   Size of grid cell in fm.
   The default is 0.2 fm, sufficient to achieve ~99.9% precision for the event properties.
   This can reasonably be increased as far as the nucleon width; beyond that and precision suffers significantly.

The grid will always be a square *N* × *N* array, with *N* = ceil(2*max/step).
So e.g. the default settings (max = 10 fm, step = 0.2 fm) imply a 100 × 100 grid.
The ceiling function ensures that the number of steps is always rounded up, so e.g. given max = 10 fm and step 0.3 fm, the grid will be 67 × 67.
In this case, the actual grid max will be marginally increased (max = nsteps*step/2).

Regardless of the collision system, the code will always approximately center the overlap region on the grid.

Configuration files
~~~~~~~~~~~~~~~~~~~
All options may be saved in configuration files and passed to the program via the ``-c, --config-file`` option.
Config files follow a simple ``key = value`` syntax, and lines beginning with a ``#`` are comments.
The key for each option is its long option without the ``--`` prefix.
Here's an example including all options::

   # specify the projectile option twice
   projectile = Pb
   projectile = Pb
   number-events = 1000

   # don't print event properties to stdout, save to HDF5
   quiet = true
   output = PbPb.hdf

   reduced-thickness = 0
   fluctuation = 1
   nucleon-width = 0.5
   cross-section = 6.4
   normalization = 1

   # leave commented out for min-bias
   # b-min =
   # b-max =

   grid-max = 10
   grid-step = 0.2

Multiple config files can be given and they will be merged, so options can be separated into modular groups.
For example, one could have a file ``common.conf`` containing settings for all collision systems and files ``PbPb.conf`` and ``pp.conf`` for specific collision systems::

   # common.conf
   reduced-thickness = 0.2
   fluctuation = 1.5
   nucleon-width = 0.6

   # PbPb.conf
   projectile = Pb
   projectile = Pb
   number-events = 10000
   grid-max = 9

   # pp.conf
   projectile = p
   projectile = p
   number-events = 100000
   grid-max = 3

To be used like so::

   trento -c common.conf -c PbPb.conf
   trento -c common.conf -c pp.conf

If an option is specified in a config file and on the command line, the command line overrides.

Examples
--------
Run a thousand lead-lead events using default settings and save the event data to file::

   trento Pb Pb 1000 > PbPb.dat

Run proton-lead events with a larger cross section (for the higher beam energy) and also compress the output::

   trento p Pb 1000 --cross-section 7.1 | gzip > pPb.dat.gz

Suppress printing to stdout and save events to HDF5::

   trento p Pb 1000 --cross-section 7.1 --quiet --output events.hdf

Uranium-uranium events at RHIC (smaller cross section) using short options::

   trento U U 1000 -x 4.2

Deformed gold-gold with an explicit nucleon width::

   trento Au2 Au2 1000 -x 4.2 -w 0.6

Simple sorting and selection (e.g. by centrality) can be achieved by combining standard Unix tools.
For example, this sorts by centrality (multiplicity) and selects the top 10%::

   trento Pb Pb 1000 | sort -rgk 4 | head -n 100

Loading data into Python
~~~~~~~~~~~~~~~~~~~~~~~~
T\ :sub:`R`\ ENTo is not designed specifically to work with Python (it is designed to be maximally flexible), but Python is extremely powerful and the authors have extensive experience using it for data analysis.

One way to load event properties is to save them to a text file and then read it with ``np.loadtxt``.
Here's a nice trick to avoid the temporary file:

.. code:: python

   import subprocess
   import numpy as np

   proc = subprocess.Popen('trento Pb Pb 1000'.split(), stdout=subprocess.PIPE)
   data = np.array([l.split() for l in proc.stdout], dtype=float)
   proc.stdout.close()

Now the ``data`` array contains the event properties.
It can be sorted and selected using numpy indexing, for example to sort by centrality as before:

.. code:: python

   data_sorted = data[data[:, 3].argsort()[::-1]]
   central = data[:100]

Text files are easily read by ``np.loadtxt``.
The header will be ignored by default, so this is all it takes to read and plot a profile:

.. code:: python

   import matplotlib.pyplot as plt

   profile = np.loadtxt('events/0.dat')
   plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)

Reading HDF5 files requires `h5py <http://www.h5py.org>`_.
Simple example:

.. code:: python

   import h5py

   h5file = h5py.File('events.hdf')
   dataset = h5file['event_0']

   # extract the grid
   profile = dataset[:]

   # read event properties
   mult = dataset.attrs['mult']
   e2 = dataset.attrs['e2']

Attribution
-----------
If you make use of this software in your research, please `cite it <http://inspirehep.net/record/1334386>`_.
The BibTeX entry is::

   @article{Moreland:2014oya,
         author         = "Moreland, J. Scott and Bernhard, Jonah E. and Bass,
                           Steffen A.",
         title          = "{An effective model for entropy deposition in high-energy
                           pp, pA, and AA collisions}",
         year           = "2014",
         eprint         = "1412.4708",
         archivePrefix  = "arXiv",
         primaryClass   = "nucl-th",
         SLACcitation   = "%%CITATION = ARXIV:1412.4708;%%",
   }

Running ``trento --bibtex`` will also print this entry.
