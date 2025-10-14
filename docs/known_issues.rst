Known Issues
------------

Using Cmake directly to rebuild the c code
==========================================
We advise against using cmake directly to rebuild the c code.
It is preferable to use ``pip install -e .[dev]`` for initial developer installation
and rebuilding the C++ code.

Using cmake directly leaves behind many cache files directories that can cause
failures when subsequently running ``pip install``. If you encounter errors like
the following when attempting to pip install lephare, remove any cmake files
except for the top level CMakeLists.txt file and try again.

.. code-block:: bash

   Error: could not load cache
        Traceback (most recent call last):
        File "/private/var/folders/j7/6btbg63d4816jknbbfxh1dfc0000gn/T/pip-build-env-f05yybu9/overlay/lib/python3.12/site-packages/setuptools/command/editable_wheel.py", line 155, in run
          self._create_wheel_file(bdist_wheel)
        File "/private/var/folders/j7/6btbg63d4816jknbbfxh1dfc0000gn/T/pip-build-env-f05yybu9/overlay/lib/python3.12/site-packages/setuptools/command/editable_wheel.py", line 357, in _create_wheel_file
          files, mapping = self._run_build_commands(dist_name, unpacked, lib, tmp)

Environmental Variables
=======================
The environment variables are set when lephare is imported by python.
If the import command is run a second time in the same session, the module is
not reinitialized, and changes to environment variables will not be propagated
through lephare. 

If you change or reset the environment variables from the command line, we
recommend restarting the kernel or python REPL to ensure that lephare correctly
reads and uses the updated environment variables.


NaNs in the output
==================
By default lephare uses -99 in place of NaNs. The ``process`` function will
handle NaNs and return them -99.

Compilers on MacOS
==================
We have seen compiler problems on macs. We hope that following the installation
instructions should avoid these issues. If you encounter segmentation faults when
using lephare on MacOS or failures to compile when running ``pip install -e .'[dev]'``
check that clang compilers are set as the default.

.. code-block:: bash

    echo $CC
    #arm64-apple-darwin20.0.0-clang
    echo $CXX
    #arm64-apple-darwin20.0.0-clang++

If these are not set to clang, you can attempt to set them using conda with the
forllowing commands:

.. code-block:: bash

    conda install -c conda-forge cxx-compiler

PyPI installation from source
=============================
When no PyPI binaries are available there may be issues installing from source.
If we have not covered your operating system or Python version pip will
attempt to install from source. If you have problems installing from source via 
PyPI it may be simpler to follow the developer installation instructions.

AVX2 compiler flags
===================
Advanced Vector eXtensions (AVX) can increase performance when they are available.
We have had an installation issue on high performance clusters where the node
that lephare was installed on had access to AVX2 and compiled with the flag set 
but then failed when it was run on a node without them due to old hardware. 
Advanced users can switch this flag off if they encounter this issue which can
be manifest in the vague error: illegal instruction (core dumped).

Installation on Fedora redhat
=============================
There may be an issue compiling lephare on Fedora/redhat, based on an attempt
on a rocky 8.8 distribution (gcc 8.5.0). The compilation error is related to the
std::filesystem standard library and is alleviated by adding the compilation flag
stdc++fs. You can see information at `stack overflow <https://stackoverflow.com/questions/71548227/undefined-reference-to-stdfilesystem-cxx11>`_.