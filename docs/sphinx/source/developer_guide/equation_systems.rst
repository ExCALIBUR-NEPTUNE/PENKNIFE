Equation Systems
=========================


The base class is the PlasmaSystem class.  This class contains the pointers to other key objects and its derived classes determine the evolution of the plasma quantities.


Magnetic Field
--------------

Every simulation must have a 3D magnetic field.
The MagneticField class handles the creation of the magnetic field from user input.
There are always three components, regardless of the dimension of the mesh - these are stored as Nektar++ expansion lists.
It also stores the unit field vectors and squared magnitude, which are used frequently throughout the code.
The magnetic field can be specified via data file or equation.
Since PENKNIFE is designed to operate in scenarios with high cylindrical symmetry (e.g tokamaks), an option is provided to specify a 2D field and have it extruded azimuthally.

The Nektar++ interpolation functionality is used to transfer field data defined on arbitrary points (perhaps imported from an eqdsk file) onto the finite element quadrature points.


Species
-------

PENKNIFE is a multispecies code, with the number and details of species being specified by the user in the config file.
These details are stored in maps in the PlasmaSystem class and can be accessed using GetIons(), GetNeutrals() or GetSpecies() for ions, neutrals and all species respectively.
Typically structured bindings are used to access key-value pairs when iterating over these maps.


Field Indexing
--------------

PENKNIFE creates Nektar++ expansion lists, hereafter referred to as "fields" for each evolved quantity.
Shared pointers to these fields are stored in an array of independent fields called m_indfields.
In order to modify or read the correct field from this array, there is an indexing system.
Each species of ion or neutral has an associated number of fields.
The index in the global array of fields can be determined from the species ID and the name of the field (e.g. n,v,e).
Generally the letters n, v and e are used in variable names to refer to (number) density, momentum or internal energy respectively.

In addition to the per-species fields, there may also be independent fields not particular to any species, for instance electron energy or vorticity.
These fields are stored after the per-species fields in m_indfields, and have their own index.
There may also be dependent fields such as electron density or momentum which are calculated from quasineutrality and current conservation considerations.
These may be stored in m_fields[0] and m_fields[1].

Field Access
------------

*The following information may change with the Nektar++ redesign*

The actual data in an expansion list/field can be accessed using GetPhys(), UpdatePhys(), GetCoeffs() and UpdateCoeffs(), depending on whether one wants read-only or write access, and physical values or coefficient values.
These functions return an Array<OneD, NekDouble>, which may be indexed using square bracket operators as with most C++ style containers.
To go from coefficient data to physical data, call the expansion list member function BwdTrans(coeffs, phys).
To go from physical data to coefficient data, call the expansion list member function FwdTrans(phys, coeffs).

Vector operations
-----------------

*The following information may change with the Nektar++ redesign*

Functions in the namespace Vmath are Nektar++ optimisations for simple operations on the Array class.
For instance Vmath::Vadd(n_pts, array_a, 1, array_b, 1, array_c, 1) adds array_a and array_b into c, where n_pts is the size of the arrays, and the 1s refer to the stride between data.
There are efforts to replace these functions where possible with fused operations that perform more complex arithmetic on a pointwise basis.
These functions will be easier to replace with GPU kernels.

Variable Conversion
-------------------

Occasionally it is desirable to transform a stored field onto another one.
One example of when this is required is to calculate temperature from internal energy and density, as is required for the diffusive part of the energy conservation equation.
For this purpose, the VariableConverter class is provided, as well as an ideal gas equation of state class.
