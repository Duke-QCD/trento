Tests & internals
=================

Unit tests
----------
.. image:: https://travis-ci.org/Duke-QCD/trento.svg?branch=master
    :target: https://travis-ci.org/Duke-QCD/trento
    :width: 90px

T\ :sub:`R`\ ENTo has a thorough `test suite <https://github.com/Duke-QCD/trento/tree/master/test>`_.
To run the tests, build the code in debug mode::

   mkdir build-debug && cd build-debug
   cmake .. -DCMAKE_BUILD_TYPE=Debug
   make tests

Tests rely on the `Catch <https://github.com/philsquared/Catch>`_ C++ test framework, which CMake downloads automatically during the build process.
If for some reason this fails, manually download the header `catch.hpp <https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp>`_ and place it in the ``test`` directory.

If `gcovr <http://gcovr.com>`_ is installed a code coverage report will be generated.

If any `sanitizers <https://github.com/google/sanitizers>`_ are available they will run.
The sanitizers check for runtime errors such as memory leaks.
They are included in recent versions of GCC and Clang.

Fitting the cross section
-------------------------
It is vital that the model reproduce the inelastic nucleon-nucleon cross section `\sigma_{NN}`.
This condition may be written

.. math::

      \sigma_{NN} = \int d^2b \, P_\text{coll}(b),

where

.. math::

   P_\text{coll}(b) = 1 - \exp[-\sigma_{gg}T_{AB}(b)]

is the probability of two nucleons colliding at impact parameter `b`.

`T_{AB}` is the overlap integral of the two nucleons' thickness functions, easily computed given Gaussian nucleons of width `w`:

.. math::

   \begin{align*}
      T_A(x, y) &= \frac{1}{2\pi w^2} \exp\biggl( -\frac{x^2 + y^2}{2w^2} \biggr), \\[1em]
      T_{AB}(b) &= \int dx \, dy \, T_A(x - b/2, y) T_B(x + b/2, y) \\
                &= \frac{1}{4\pi w^2} \exp\biggl( -\frac{b^2}{4w^2} \biggr).
   \end{align*}

`\sigma_{gg}` is an effective parton-parton cross section now determined by the relation

.. math::

   \sigma_{NN} = \int_0^{b_\text{max}} 2\pi b \, db \,
      \biggl\{
         1 - \exp\biggl[
            -\frac{\sigma_{gg}}{4\pi w^2}
            \exp\biggl( -\frac{b^2}{4w^2} \biggr)
         \biggr]
      \biggr\},

where `b_\text{max}` is the maximum impact parameter for a collision.
Let `b_\text{max} = Aw`, i.e. some number of nucleon widths (the actual code uses `A = 6`).

After appropriate change of variables, this relation may be written

.. math::

   \frac{\sigma_{NN}}{4\pi w^2} =
      \frac{A^2}{4} +
      \text{Ei}\biggl( -e^{-A^2/4} \frac{\sigma_{gg}}{4\pi w^2} \biggr) -
      \text{Ei}\biggl( -\frac{\sigma_{gg}}{4\pi w^2} \biggr)

where `\text{Ei}` is the `exponential integral <https://en.wikipedia.org/wiki/Exponential_integral>`_.
This is still a transcendental equation but it can be quickly solved numerically for a given cross section and nucleon width.

The code determines `\sigma_{gg}` at runtime by solving this equation using a standard root finding algorithm.
The relevant function is ``compute_cross_sec_param`` in `nucleon.cxx <https://github.com/Duke-QCD/trento/blob/master/src/nucleon.cxx>`_.
Note that it actually optimizes over the dimensionless variable `\log({\sigma_{gg}/4\pi w^2})`.

The nucleon unit test verifies that the cross section is accurately reproduced.

.. note::

   The semi-analytic cross section procedure described above only works for a collision involving two Gaussian particles.
   If nucleon substructure is enabled, i.e. constituent number *m > 1*, the code uses a Monte Carlo algorithm to fit the cross section parameter `\sigma_{gg}` numerically.
   This procedure is somewhat slow, and it scales poorly with the number of constituents.
   For *m > 10*, it could easily take a minute or longer, so the code caches the value of `\sigma_{gg}` in the user's XDG_DATA_HOME directory.

List of classes
---------------
This section is automatically generated from the source code by `Doxygen <http://www.stack.nl/~dimitri/doxygen>`_ and converted into `Sphinx <http://sphinx-doc.org>`_ format by `Breathe <https://breathe.readthedocs.org>`_.
In order to emphasize functionality---rather than implementation details---private class methods are not shown.

.. highlight:: cpp

Collider
~~~~~~~~
.. doxygenclass:: trento::Collider
   :members: Collider, run_events

Event
~~~~~
.. doxygenclass:: trento::Event
   :members:

Output
~~~~~~
.. doxygenclass:: trento::Output
   :members:

Nucleus
~~~~~~~
.. doxygenfunction:: trento::Nucleus::create

.. doxygenclass:: trento::Nucleus
   :members: radius, sample_nucleons, ~Nucleus
   :protected-members:

Nucleus types
'''''''''''''
.. doxygenclass:: trento::Proton
.. doxygenclass:: trento::Deuteron
.. doxygenclass:: trento::WoodsSaxonNucleus
.. doxygenclass:: trento::DeformedWoodsSaxonNucleus
.. doxygenclass:: trento::ManualNucleus

NucleonData
~~~~~~~~~~~
.. doxygenclass:: trento::NucleonData
   :members: NucleonData, is_participant, constituents_exist, x, y, z, set_position

NucleonCommon
~~~~~~~~~~~~~~~
.. doxygenclass:: trento::NucleonCommon
   :members: NucleonCommon, max_impact, boundary, thickness, participate, sample_constituent_positions, set_participant

Fast exponential
~~~~~~~~~~~~~~~~
.. doxygenclass:: trento::FastExp
   :members:
