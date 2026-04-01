====================================
Braginskii Closure
====================================

Collision Frequencies
---------------------

The first step in calculating the Braginskii heat fluxes and frictions is to calculate collsiion frequencies between all pairs of species.

The equation for collision frequency of species a with species b is:

.. math:: \nu_{ab} = \frac{Z_a^2 Z_b^2 e^4 \ln{\Lambda_{ab}}}{3\pi^{3/2}\epsilon_0^2} \frac{n_b}{\left(v_a^2+v_b^2\right)^{3/2}} \frac{1}{\mu_{ab}m_a}

Note that the collision rate is for Coulomb scattering events, thus :math:`\nu_{ab}\neq\nu_{ba}`.
:math:`nu_{ab}` is the frequency at which species a undergoes Coulomb scattering due to collisions with species b.

In fact, due to conservation of momentum:

.. math:: m_a n_a \nu_{ab} = m_b n_b \nu_{ba}


Coulomb Logarithm
-----------------

The Coulomb logarithm :math:`\ln\Lambda` appears in the collision frequency equation in order to account for the cumulative effect of multiple small-angle collisions.
Its value depends on the species involved in the collisions, and has qualitatively different forms for electron and ion collisions.

.. math:: \ln{\Lambda_{ee}} =   30.4 - \frac{1}{2}\ln{n_e} + \frac{5}{4}\ln{T_e} - \sqrt{10^{-5}+\left(\ln{T_e}-2\right)^2/16}
.. math:: \ln{\Lambda_{ei}} =   \begin{cases}
                                    10                                                                      & \text{if}\ T_e         \lt0.1eV\       \text{or}\ n_e \lt 10^{10}m^{-3} \\
                                    30 - \frac{1}{2}\ln{n_e} - \ln{Z} + \frac{3}{2}\ln{T_e}                 & \text{if}\ T_i m_e/m_i \lt T_e         \lt 10Z^2                        \\
                                    31 - \frac{1}{2}\ln{n_e} + \ln{T_e}                                     & \text{if}\ T_i m_e/m_i \lt 10Z^2       \lt T_e                          \\
                                    23 - \frac{1}{2}\ln{n_e} + \frac{3}{2}\ln{T_e} - \ln{\left(Z^2A\right)} & \text{if}\ T_e         \lt T_i m_e/m_i
                                \end{cases}
.. math:: \ln{\Lambda_{ii}} =   29.91 - \ln{\left[\frac{Z_aZ_b\left(m_a+m_b\right)}{m_aT_b+m_bT_a} \sqrt{\frac{n_aZ_a^2}{T_a} + \frac{n_bZ_b^2}{T_b}}\right]}
        

Heat flux
---------

The collisional parallel heat flux for species i is given by:

.. math:: \mathbf{q}_i = C_i\frac{n T_i}{m_i \nu_i} \left(\mathbf{b}\cdot\nabla T_i\right)\mathbf{b}

where

.. math:: \nu_i = \sum_j \nu_{ij}

:math:`C_{ion} = 3.9` and  :math:`C_{electron} = 3.19`

Friction
--------

Momentum is exchanged between two species by collisions according to:

.. math:: F_{ab} = \nu_{ab} m_a n_a \left(v_b - v_a\right)

Associated with this frictional momentum exchange is a frictional heating (Joule heating):

.. math:: Q_{ab} = \frac{m_b}{m_a + m_b} \left(v_b - v_a\right) F_{ab}

Both species' internal energy increases as a result of collisions.

Heat Exchange
-------------

Two collisional species exchange heat if they are at different temperatures

The heat exchanged from species b to species a by this mechanism is:

.. math:: W_{ab} = 3 \nu_{ab} n_a \frac{m_a}{m_a+m_b} \left(T_b - T_a\right)