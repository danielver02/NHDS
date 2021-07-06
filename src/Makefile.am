if WITH_HDF5
AM_FCFLAGS = -I$(HDF5_INC) -fbounds-check -ffpe-trap=invalid,zero,overflow -O2 -ffixed-line-length-none -ffree-line-length-none -Wunused -lm -ldl
NHDS_LDADD = -L$(HDF5_LIB) $(HDF5_LDFLAGS)
else
AM_FCFLAGS = -fbounds-check -ffpe-trap=invalid,zero,overflow -O2 -ffixed-line-length-none -ffree-line-length-none -Wunused -lm -ldl
endif


bin_PROGRAMS = NHDS

if WITH_HDF5
NHDS_SOURCES = NHDS_with_hdf5.f90 write_delta_f.f90
else
NHDS_SOURCES = NHDS.f90
endif

NHDS_SOURCES += calc_polarization.f90 calc_xi.f90 dispfunct.f90 disptensor.f90 newton_method.f90 parameters.f90 suscep.f90 waveenergy.f90
