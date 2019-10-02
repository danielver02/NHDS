AM_FCFLAGS = -fbounds-check -O2 -ffixed-line-length-none -ffree-line-length-none -Wunused -lm -I$(HDF5_INC)
AM_LDFLAGS = -ldl -L$(HDF5_LIB) $(HDF5_LDFLAGS)

bin_PROGRAMS = NHDS
NHDS_SOURCES = NHDS.f90 calc_polarization.f90 calc_xi.f90 dispfunct.f90 disptensor.f90 newton_method.f90 parameters.f90 suscep.f90 waveenergy.f90 write_delta_f.f90
