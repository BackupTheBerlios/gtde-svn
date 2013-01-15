# First check if we have the fortrand code
if [ ! -f m1qn3.f ]
then
	echo "The file m1qn3.f does not exist. Please download it at:"
	echo "https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html"
	echo "!!!!    Error!!!!"
	exit
fi
# Compile the fortrand code in a shared library
echo "compiling fortrand"
gfortran -shared -fPIC m1qn3.f -o libm1qn3.so
echo "-- Fortran code into shared library."
# Patch m1qn3mex.c
if [ ! -f m1qn3mex.c ]
then
	echo "The file m1qn3mex.c does not exist. Please download OptiTool at:"
	echo "http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php"
	echo "And find that file in the Solvers/m1qn3/Source folder."
	echo "!!!!    Error!!!!"
	exit
fi
# Patch the c file
patch m1qn3mex.c m1qn3mex.diff
echo "-- File m1qn3mex.c patched."
# Compile wrapper
echo "compiling wrapper with mex"
mex CFLAGS='$CFLAGS -std=c99' m1qn3mex.c -L. -lm1qn3
# Final message
echo "-- Wrapper compiled."
echo "-- Please do not forgot to add $(pwd) to the LD_LIBRARY_PATH before starting MATLAB."
