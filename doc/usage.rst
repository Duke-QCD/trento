Usage
=====

.. highlight:: shell

Usage::

   trento [options] projectile projectile [number-events = 1]

Projectiles must be standard symbols, e.g. ``Pb`` for lead and ``p`` for proton.
This would run one thousand proton-lead events with default settings::

   trento p Pb 1000

The default output mode is to print these event properties to stdout::

   event_number total_entropy impact_param Npart e_2 e_3 e_4

where the ``e_n`` are eccentricity harmonics.
Give the ``-q/--quiet`` option to disable this output.

Event profiles are not output by default; the ``-o/--output`` option must be explicitly given.
If an HDF5 filename is given (``.h5``, ``.hd5``, ``.hdf5``), all profiles will be written to file.
Otherwise, the option is interpreted as a directory and profiles will be written to numbered text files in the directory.

All options (including projectiles and number of events) can be saved in configuration files and used via the ``-c/--config-file`` option.
If an option is specified in both a config file and on the command line, the command line overrides the config file.
