======================
Particle Specification
======================

PENKNIFE can use the NESO-Particles and NESO frameworks to allow users to use computational particles which interact with the fluid solver.

A Particle system is enabled using the :code:`<PARTICLES>` field in the XML config file.

Global parameters for the particle system are specified in :code:`<PARAMETERS>`

Particle information is specified on a per-species basis, within :code:`<SPECIES>`

Each species must have mass and charge fields, as well as initial conditions, and may optionally have sources and sinks.

Initial Conditions
------------------

The species initial conditions are specified in the :code:`<INITIAL N="num">` field.

The density field can be specified by an analytic function for a probability distribution (everywhere between 0 and 1).
During initialisation, `num` particles are created using adaptive rejection sampling over this distribution.

Sources
-------

Sources are specified under a :code:`<SOURCE N="num">` field.  Multiple sources can be created for each species.
Here, `num` indicates the number of particles added at each timestep.

Sinks
-----

Sinks are specified under a :code:`<SINK>` field.  Multiple sinks can be created for each species.
The specified function indicates the probability that a particle is removed at each timestep.
The probability function is evaluated at each particle's position, which is then marked for removal with that probability.

Boundary Conditions
-------------------

The particles may interact with the boundary regions by either reflecting, being absorbed, or taking part in a surface reaction.



