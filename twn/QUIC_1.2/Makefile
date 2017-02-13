# Author: Matyas A. Sustik
#
# Compilation for QUIC written by Cho-Jui Hsieh and Matyas A. Sustik based
# on:
#
#   Sparse Inverse Covariance Matrix Estimation Using Quadratic
#   Approximation. Cho-Jui Hsieh, Matyas A. Sustik, Inderjit S. Dhillon,
#   Pradeep Ravikumar, Advances in Neural Information Processing Systems,
#   vol. 24, 2011, p. 2330--2338.

# Question? Bug report? Please send email to msustik@gmail.com.

# You may compile QUIC using the MEX compiler provided by Matlab,
# and that is the recommended approach for Windows and MacOS users. See
# the README file for details.

# Please compile the Octave executable by issuing at a shell (not
# within Octave):
#
# > make QUIC.mex
#
# To compile the Matlab executable, first specify the location of your
# Matlab installation directory by setting the MATLAB_ROOT environment
# variable (below), then execute in a shell (not from within Matlab):
#
# > make QUIC.mexa64
#
# which produces a 64-bit MEX file. To compile for 32-bit Matlab on
# Linux, please use:
#
# > make QUIC.mexglx
#
# On MacOS and Windows use the appropriate extension.
#
#
# DEBUGGING
#
# You may compile a debug version of the program which takes care of
# attaching gdb and popping up a debug window, without the need to modify
# your Matlab or Octave environment.
#
# For example, when working under Octave, compile the debug version with:
#
# > make QUIC_g.mex
#

# A gdb debug window should pop up when QUIC_g is called. You will
# need startgdb.h that we make available separately and upon request.
#
# Please report the bug and the fix if you have found one.

# Contact your administrator if you do not know where your matlab is
# installed.
MATLAB_ROOT = /lusr/share/software/matlab-2011a

# Normally you would not need to modify lines bellow this one.  If you
# need to, please let us know how we can improve this Makefile.

MATLAB_INCLUDE = -I$(MATLAB_ROOT)/extern/include
# Octave has a nice feature to tell us the library and include locations:
OCTAVE_INCLUDE = $(shell mkoctfile -p INCFLAGS)
OCTAVE_LIBS = $(shell mkoctfile -p LFLAGS)

CXX=g++

CXXFLAGS = -Wall -pedantic -ansi -fpic -pthread -shared  -fno-omit-frame-pointer -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64

CXXOPTFLAGS = -O2 -DNDEBUG $(CXXFLAGS)
CXXDBGFLAGS = -g -DNDEBUG -DGDEBUG $(CXXFLAGS)

LDMEXFLAGS = -Wall -pthread -shared -lmwblas -lmwlapack -lmx -lmex -lmat -lm -lgfortran

# Linker flags for Matlab on Linux:
LDFLAGS_32 = -m32 -Wl,--version-script,$(MATLAB_ROOT)/extern/lib/glnx86/mexFunction.map -Wl,-rpath-link,$(MATLAB_ROOT)/bin/glnx86 -L$(MATLAB_ROOT)/bin/glnx86 $(LDMEXFLAGS)
LDFLAGS_64 = -m64 -Wl,--version-script,$(MATLAB_ROOT)/extern/lib/glnxa64/mexFunction.map -Wl,-rpath-link,$(MATLAB_ROOT)/bin/glnxa64 -L$(MATLAB_ROOT)/bin/glnxa64 $(LDMEXFLAGS)

# Linker flags for Matlab on MacOS:
LDMACFLAGS_32 = -m32 -Wl,-exported_symbols_list,$(MATLAB_ROOT)/extern/lib/maci32/mexFunction.map -Wl,-rpath,$(MATLAB_ROOT)/bin/maci32 -L$(MATLAB_ROOT)/bin/maci32 $(LDMEXFLAGS)
LDMACFLAGS_64 = -m64 -Wl,-exported_symbols_list,$(MATLAB_ROOT)/extern/lib/maci64/mexFunction.map -Wl,-rpath,$(MATLAB_ROOT)/bin/maci64 -L$(MATLAB_ROOT)/bin/maci64 $(LDMEXFLAGS)

# Linker flags for Octave on Linux:
LDOCTMEXFLAGS = -shared -Wl,-Bsymbolic -loctave -llapack -lblas -lfftw3 -lfftw3f -lreadline -lncurses -ldl -lz -lm -lgfortranbegin -lgfortran $(OCTAVE_LIBS)

.SUFFIXES:

# C++ rules for MEX wrappers:
%-mex_32_g.o: %-mex.cpp
	$(CXX) $(CXXDBGFLAGS) $(MATLAB_INCLUDE) -m32 -c $< -o $@

%-mex_32.o: %-mex.cpp
	$(CXX) $(CXXOPTFLAGS) $(MATLAB_INCLUDE) -m32 -c $< -o $@

%-mex_64_g.o: %-mex.cpp
	$(CXX) $(CXXDBGFLAGS) $(MATLAB_INCLUDE) -m64 -c $< -o $@

%-mex_64.o: %-mex.cpp
	$(CXX) $(CXXOPTFLAGS) $(MATLAB_INCLUDE) -m64 -c $< -o $@

%-oct_g.o: %-mex.cpp
	$(CXX) $(CXXDBGFLAGS) $(OCTAVE_INCLUDE) -c $< -o $@

%-oct.o: %-mex.cpp
	$(CXX) $(CXXOPTFLAGS) $(OCTAVE_INCLUDE) -c $< -o $@

# C++ rules for algorithm programs:
%_32_g.o: %.cpp
	$(CXX) $(CXXDBGFLAGS) -m32 -c $< -o $@

%_32.o: %.cpp
	$(CXX) $(CXXOPTFLAGS) -m32 -c $< -o $@

%_64_g.o: %.cpp
	$(CXX) $(CXXDBGFLAGS) -m64 -c $< -o $@

%_64.o: %.cpp
	$(CXX) $(CXXOPTFLAGS) -m64 -c $< -o $@

%_g.o: %.cpp
	$(CXX) $(CXXDBGFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXOPTFLAGS) -c $< -o $@

%_p.o: %.cpp
	$(CXX) $(CXXOPTFLAGS) -ggdb -c $< -o $@

# Link Matlab executables:
%_g.mexa64 : %-mex_64_g.o %_64_g.o
	$(CXX) $(LDFLAGS_64) $^ -o $@

%.mexa64 : %-mex_64.o %_64.o
	$(CXX) $(LDFLAGS_64) $^ -o $@

%_g.mexglx : %-mex_32_g.o %_32_g.o
	$(CXX) $(LDFLAGS_32) $^ -o $@

%.mexglx : %-mex_32.o %_32.o
	$(CXX) $(LDFLAGS_32) $^ -o $@

%_g.mexmaci64 : %-mex_64_g.o %_64_g.o
	$(CXX) $(LDMACFLAGS_64) $^ -o $@

%.mexmaci64 : %-mex_64.o %_64.o
	$(CXX) $(LDMACFLAGS_64) $^ -o $@

%_g.mexmaci32 : %-mex_32_g.o %_32_g.o
	$(CXX) $(LDMACFLAGS_32) $^ -o $@

%.mexmaci32 : %-mex_32.o %_32.o
	$(CXX) $(LDMACFLAGS_32) $^ -o $@

# Link Octave executables:
%_g.mex : %-oct_g.o %_g.o
	g++ $(LDOCTMEXFLAGS) $^ -o $@

%.mex : %-oct.o %.o
	g++ $(LDOCTMEXFLAGS) $^ -o $@ 

.SECONDARY :

clean :
	rm -f *.o *.mex* *.oct


