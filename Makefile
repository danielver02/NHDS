#This file is part of NHDS.
#
#    NHDS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NHDS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NHDS.  If not, see <http://www.gnu.org/licenses/>.

SRCS := $(shell ls src/*.f90)
OBJS := $(SRCS:src/%.f90=obj/%.o)
HDF5 := /usr/local/hdf5
FORTRANLIB := -I$(HDF5)/include $(HDF5)/lib/libhdf5_fortran.a
LIBSHDF := $(FORTRANLIB) $(HDF5)/lib/libhdf5.a
LIBZ :=/usr/local/lib/libz.a

COMPILER := gfortran -fbounds-check -O2 -ffixed-line-length-none -ffree-line-length-none -Wunused -lm


bin/NHDS: $(OBJS) 
	$(COMPILER)  $(OBJS) -o $@ $(LIB64)  $(LIBSHDF) $(LIBZ)


obj/%.o:src/%.f90
	$(COMPILER)  -c  $< -o $@ 


clean:
	rm -f $(OBJS) bin/NHDS

