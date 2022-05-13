====Dependencies ====

  - scons (http://www.scons.org/)

Compile-time libraries for precomputed table of Fresnel Integrals:

  - MPFR (http://www.mpfr.org/)
  - GMP (https://gmplib.org/)

For test programs (will be ignored if missing):

  - cairo (http://cairographics.org/)
  - vxl (http://vxl.sourceforge.net/)

All can be install on Debian/Ubuntu with:

  sudo apt-get install scons libgmp-dev libmpfr-dev libcairo2-dev libvxl1-dev

==== Building ====

libcornu can be build by running the `scons` command in the root directory of
this distribution.  The computation of the Fresnel table may take some time,
e.g. 71 minutes on 7 threads on a Intel(R) Core(TM) i7-2600 CPU @ 3.40GHz.

==== Running ====

If cairo is available, a test program will be built to 'bin/test-cornu'. It
optionally takes an integer as argument to be used as a random seed.  The
following calls are valid:

./bin/test-cornu
./bin/test-cornu 24
./bin/test-cornu time

The last call uses the current UNIX time as the integer.  It will plot the
resulting cornu spiral completion and the input constraints to 'cornu.pdf'
in the current directory.

==== Linking ====

Static and shared libraries will be built in the 'lib/' directory.

Necessary include files will be placed in 'include/vptree/'
