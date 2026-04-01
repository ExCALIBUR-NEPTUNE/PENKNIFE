User Input
=========================

.. role:: xml(code)
   :language: xml

The details of each simulation are specified in two xml files: the mesh and the config.



Particle Specification
======================

PENKNIFE can use the NESO-Particles and NESO frameworks to allow users to use computational particles which interact with the fluid solver.

A Particle system is enabled using the :xml:`<PARTICLES>` field in the XML config file.

Global parameters for the particle system are specified in :xml:`<PARAMETERS>`

Particle information is specified on a per-species basis, within :xml:`<SPECIES>`

Each species must have mass and charge fields, as well as initial conditions, and may optionally have sources and sinks.

Initial Conditions
------------------

The species initial conditions are specified in the :xml:`<INITIAL N="num">` field.

The density field can be specified by an analytic function for a probability distribution (everywhere between 0 and 1).
During initialisation, `num` particles are created using adaptive rejection sampling over this distribution.

Sources
-------

Sources are specified under a :xml:`<SOURCE N="num">` field.  Multiple sources can be created for each species.
Here, `num` indicates the number of particles added at each timestep.

The following xml indicates that there is a point source of particles at (2.75, 0), adding 100 particles of weight 0.05 at each time step, with thermal velocities at temperature 1eV

.. code-block:: xml

    <SOURCES>
        <P N="100">
            <E VAR="W" VALUE="0.05"/>
            <E VAR="X" VALUE="2.75"/>
            <E VAR="Y" VALUE="0"/>
            <E VAR="T" VALUE="1"/>
        </P>
    </SOURCES>

Sinks
-----

Sinks are specified under a :xml:`<SINK>` field.  Multiple sinks can be created for each species.
The specified function indicates the probability that a particle is removed at each timestep.
The probability function is evaluated at each particle's position, which is then marked for removal with that probability.

Boundary Conditions
-------------------

The particles may interact with the boundary regions by either reflecting, being absorbed, or taking part in a surface reaction.

Optionally, the user can use the Reactions/VANTAGE framework by populating the :xml:`<VANTAGE>` XML field.

The three currently supported reactions are ionisation, recombination and charge exchange.


Reactions coupling
==================

Ionisation
----------

Recombination
-------------

Charge Exchange
---------------
