Class list
==========

Collider
--------
.. doxygenclass:: trento::Collider
   :members: Collider, run_events

Event
-----
.. doxygenclass:: trento::Event

Output
------
.. doxygenclass:: trento::Output

Nucleus
-------
.. doxygenfunction:: trento::Nucleus::create

.. doxygenclass:: trento::Nucleus
   :members: radius, sample_nucleons, ~Nucleus
   :protected-members:

Nucleus types
~~~~~~~~~~~~~

Proton
''''''
.. doxygenclass:: trento::Proton

Woods-Saxon
'''''''''''
.. doxygenclass:: trento::WoodsSaxonNucleus

Deformed Woods-Saxon
''''''''''''''''''''
.. doxygenclass:: trento::DeformedWoodsSaxonNucleus

Nucleon
-------
.. doxygenclass:: trento::Nucleon
   :members: Nucleon, x, y, is_participant, set_position, set_participant

Nucleon profile
---------------
.. doxygenclass:: trento::NucleonProfile

Fast exponential
----------------
.. doxygenclass:: trento::FastExp
