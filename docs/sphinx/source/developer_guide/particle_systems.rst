Particle Systems
=========================

Particle Loops
--------------

Particle loops are a NESO-Particles feature that allow iteration over particle groups and subgroups.
In PENKNIFE they are used explicitly in a few places, most importantly in the inner integration routine to push or advect the particles to their next position.


Projection and Evaluation
-------------------------

Projection is the process of converting particle data into a finite element representation.
A particular particle quantity is first considered to be a function of a series of delta functions at the particle positions.
The finite element representation is then given by the polynomial coefficients which best reproduce this function.

Evaluation is conceptually simpler - it is the evaluation of the physical field values at the particle positions.

Projection is used to transfer source terms from particles to fields, as well as producing diagnostic fields for particle density.
Evaluations are used wherever fluid quantities must be used by the particles, for instance in calculating reaction rates.


VANTAGE
-------

VANTAGE capability is provided in the ReactionSystem class.
Here the reactions specified by the user under the <VANTAGE> xml tag are constructed and added to the appropriate reaction_controller.
There is one reaction controller for each reaction that occurs in the bulk, as well as one reaction controller for each boundary region identified in the xml.
The boundary reaction controllers operate in surface mode.

The reactions themselves are built from reaction kernels and reaction data.
Support functions for building reactions is given in Reactions.hpp.
Reactions may need to be templated on the mesh dimension and/or the number of velocity components.
For instance, the specular reflection reaction in 2D cylindrical coordinates is a different instantiation of the template to 2D Cartesian coordinates, since there are 3 velocity dimensions in this case.


Boundaries
----------

At boundaries the particles undergo a transformation such as specular reflection or absorption.
It is necessary to adjust the velocity of a particle that is headed for a boundary, so that it does not leave the domain.
The weight of the particle may also be decreased due to absorption
The boundary handling process goes as follows:

1.  All particles that are on a course to intersect a boundary region are identified and their present positions stored.
2.  Particles are advected until the end of the timestep, or until they hit a boundary, whichever comes first.
3.  The remaining timestep is stored for colliding particles
4.  Those which hit a boundary have the velocity and weight adjustments made.
5.  If new particles were created by a collision they are added to the particle group.

It is possible that this process must occur multiple times within one particle timestep if the particle hits a corner.
A while loop containing these steps runs until there are no particles on course to hit a boundary.