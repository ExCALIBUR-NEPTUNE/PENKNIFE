====================================
Double diffusive field with neutrals
====================================

Solves a system of :math:`2+2S` fields representing :math:`S` ion species.
This equation system is used by specifying

.. code-block:: xml

    <I PROPERTY="EQTYPE" VALUE="DoubleDiffusiveField"/>

in the :code:`<SOLVERINFO>` field.

Here the density field from the single field system is used alongside an energy field,
from which is derived the temperature.

We will ignore all advective and field terms for the first pass, and
write a pressure equation as

.. math:: \frac{3}{2} \frac{\partial p}{\partial t} + \nabla \cdot \vec{q} = S_E + Q,

where :math:`p = nk_BT` and the flux:

.. math:: \vec{q} = \kappa \cdot \nabla T

with the conductivity tensor given as:

.. math:: \kappa =  \kappa_\parallel(T) \vec{b}\vec{b} + \kappa_\perp(T) (I-\vec{b}\vec{b})

with parallel and perpendicular conductivities being functions of
temperature.

Energy sinks :math:`S_E` should again be supplied by VANTAGE, with
some heating :math:`Q` in the core being an input parameter.

Electrons are always present, and at least one ion species must be present.
Explicitly the fields are:

.. math:: \left[n_e,  \frac{3}{2}p_e, n_i, \frac{3}{2}p_i, ...\right]

The continuity equations are:

.. math:: \frac{\partial n}{\partial t} + \nabla \cdot \vec{\Gamma} = S_n
.. math:: \frac{3}{2} \frac{\partial p}{\partial t} + \nabla \cdot \vec{q} = S_E + Q,


Boundary Conditions
-------------------

The boundary conditions for this system of equations are as follows.

The energy outflow is set to:

.. math:: \vec{q}\cdot \vec{b} = \gamma nkTc_s,

The density outflow is the same as for the single field:
.. math:: \vec{\Gamma}\cdot \vec{b} = n c_s

The two options for :math:`c_s` are:
All ion species enter the sheath at the system sound velocity.

.. math:: c_s = \sqrt{\sum_i \frac{n_{0i} k_B T_e}{n_{0e} m_i}}

The ion species enter the sheath at their particular sound velocity:

.. math:: c_{si} = \sqrt{\frac{n_{0i} k_B T_e}{n_{0e} m_i}}
