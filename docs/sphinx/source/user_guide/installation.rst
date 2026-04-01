Installation
=========================

The easiest way to build the PENKNIFE is to use the `Spack package manager <https://spack.readthedocs.io/en/latest/index.html>`_.
The build process has been tested with **Spack v1.2.0.dev0**; later versions may also work, earlier than **v1.0.0** is not recommended.

Get Spack
---------
Install Spack via the `official instructions <https://spack.readthedocs.io/en/latest/getting_started.html#installation>`_.
On Debian-based systems (e.g. Ubuntu) the following should work:

.. code-block::  bash

    git clone --depth=2 https://github.com/spack/spack.git

Initialise Spack
~~~~~~~~~~~~~~~~

.. code-block::  bash

    source $SPACK_ROOT/share/spack/setup-env.sh
    # Optionally modify .bashrc such that spack is always initialised in new shells (this only has to be done once)
    echo 'source $SPACK_ROOT/share/spack/setup-env.sh' >> $HOME/.bashrc

Install intel compilers
-----------------------
*Only required for the oneapi build*

.. code-block::  bash

    spack install intel-oneapi-compilers
    spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/bin


Or if they are already installed elsewhere and version < 2024, add them

.. code-block::  bash

    spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/linux/bin/intel64
    spack compiler add `spack location -i intel-oneapi-compilers`/compiler/latest/linux/bin

N.B. There is a known problem/conflict with the mkl libraries if the compilers come from the oneapi toolkit - it is recommended not to use this

Get PENKNIFE
------------

.. code-block::  bash

    cd ..
    git clone https://github.com/ExCALIBUR-NEPTUNE/PENKNIFE.git
    cd PENKNIFE
    git submodule update --init --recursive
    source ./scripts/activate.sh
    spack install

N.B.

.. code-block::  bash

    spack install -j <nproc>

can be used to speed up compilation

Note that some of the dependencies (particularly nektar++) can take some time to install and have high memory usage.

To leave the spack environment at any point run:

.. code-block::  bash

    deactivate

