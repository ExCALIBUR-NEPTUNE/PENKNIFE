Post-Processing and Visualisation
=================================

The Nektar++ FieldConvert module is used to convert the .chk files produced by the run into .vtu files that can be viewed in Paraview.


.. code-block:: bash

    # Convert to vtu using FieldConvert (requires xml files as args)
    cd runs/[EquationSystem]/[Mesh]/[Dimension]/[Example]

    FieldConvert [config_xml] [mesh_xml] [chk_name] [vtu_name]
    # e.g. FieldConvert single_field.xml mastu.xml single_field_100.chk single_field.vtu
    paraview [vtu_name]

To convert multiple .chk files at once:

.. code-block:: bash

    for i in {0..1000}; do FieldConvert -f s[config_xml] [mesh_xml] [chk_name]_$i [vtu_name]_$i.vtu; done;
    # e.g. for i in {0..1000}; do FieldConvert -f single_field.xml mastu.xml single_field_$i.chk single_field_$i.vtu; done;

Particles can also be viewed in Paraview by opening "particle_trajectopry.h5part"
