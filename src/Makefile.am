NHDS_LDADD = -lm -ldl

if WITH_HDF5
AM_FCFLAGS = -I$(HDF5_INC) -O4 -fbounds-check -ffpe-trap=invalid,zero,overflow -ffixed-line-length-none -ffree-line-length-none -Wunused  -g -cpp -DWITH_HDF5
NHDS_LDADD += -L$(HDF5_LIB) $(HDF5_LDFLAGS)
else
AM_FCFLAGS = -fbounds-check -O4 -ffpe-trap=invalid,zero,overflow -ffixed-line-length-none -ffree-line-length-none -Wunused -g -cpp
endif

bin_PROGRAMS = NHDS

if WITH_HDF5
NHDS_SOURCES = NHDS.f90 write_delta_f.f90
else
NHDS_SOURCES = NHDS.f90
endif

NHDS_SOURCES += calc_polarization.f90 calc_xi.f90 dispfunct.f90 disptensor.f90 newton_method.f90 parameters.f90 suscep.f90 waveenergy.f90
