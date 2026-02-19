====================================
Single diffusive field with neutrals
====================================

Solves a system of :math:`1+S` fields representing :math:`S` ion species.
This equation system is used by specifying

.. code-block:: xml

    <I PROPERTY="EQTYPE" VALUE="SingleDiffusiveField"/>

in the :code:`<SOLVERINFO>` field

Electrons are always present, and at least one ion species must be present.
Explicitly the fields are:

.. math:: \left[n_e, n_i, ...\right]

The continuity equation is:

.. math:: \frac{\partial n}{\partial t} + \nabla \cdot \vec{\Gamma} = S_n

with the diffusive flux given by:

.. math::`\vec{\Gamma} = D \cdot \nabla n`

with D being the anisotropic diffusion tensor:

.. math:: D = \vec{b}\vec{b} k_\parallel + (I - \vec{b}\vec{b}) k_\perp

where the :math:`k_\parallel` and :math:`k_\perp` have the option
of being functions of fields, and :math:`\vec{b}` is the local unit
vector parallel to the magnetic field.

The source :math:`S_n` is given by VANTAGE.
The system is assumed to be isothermal, with all species at the same temperature.


Boundary Conditions
-------------------

The available boundary conditions for this system of equations are as follows.

System Bohm

All ion species enter the sheath at the system sound velocity :math:`c_s`.

.. math:: c_s = \sqrt{\sum_i \frac{n_{0i} k_B T_e}{n_{0e} m_i}}

The boundary condition should be set to Neumann in the directions
perpendicular to the magnetic field, while the parallel component of the
flux at the boundary

.. math:: \vec{\Gamma}\cdot \vec{b} = n c_s

where :math:`c_s` in this case is a constant (isothermal Bohm speed
:math:`c_s = \sqrt{kT/m_i}`).

.. math:: \frac{\partial n}{\partial \mathbf{n}} = \mp \frac{n}{c_s} \frac{\partial v_{\parallel i}}{\partial \mathbf{n}}

Species Bohm

Dirichlet

Neumann