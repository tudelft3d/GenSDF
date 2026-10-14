# FC/GPU are overridable, e.g.:
#   make GPU=1 FC=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/bin/mpif90
FC        = mpif90
GPU       ?= 0
FFLAGS_DBG= -O0 -g -fbounds-check -ffree-line-length-none -cpp	#Bounds check only in debug mode
LDFLAGS   =
EXE     = gensdf

ifeq ($(GPU),1)
    # nvfortran: no -ffree-line-length-none (gfortran-only flag, nvfortran has no line-length cap)
    FFLAGS  = -O3 -cpp -acc -gpu=cc86,rdc -Minfo=accel -DGPU_SDF
    LDFLAGS = -acc -gpu=cc86,rdc
else
    FFLAGS = -O3 -ffree-line-length-none -cpp
endif

# Compilation order respects Fortran module dependencies:
#   utils_io      -> flood_fill, spatial_hash, fast_sweep
#   spatial_hash  -> narrowband, raycast
#   all modules   -> main
SRCS = utils_io.f90 flood_fill.f90 spatial_hash.f90 narrowband.f90 raycast.f90 fast_sweep.f90 main.f90
OBJS = $(SRCS:.f90=.o)
MODS = utils_io.mod flood_fill_mod.mod spatial_hash_mod.mod \
       narrowband_mod.mod raycast_mod.mod fast_sweep_mod.mod

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) -o $@

utils_io.o: utils_io.f90
	$(FC) $(FFLAGS) -c $< -o $@

flood_fill.o: flood_fill.f90 utils_io.o
	$(FC) $(FFLAGS) -c $< -o $@

spatial_hash.o: spatial_hash.f90 utils_io.o
	$(FC) $(FFLAGS) -c $< -o $@

narrowband.o: narrowband.f90 utils_io.o spatial_hash.o
	$(FC) $(FFLAGS) -c $< -o $@

raycast.o: raycast.f90 utils_io.o spatial_hash.o
	$(FC) $(FFLAGS) -c $< -o $@

fast_sweep.o: fast_sweep.f90 utils_io.o
	$(FC) $(FFLAGS) -c $< -o $@

main.o: main.f90 utils_io.o flood_fill.o spatial_hash.o narrowband.o raycast.o fast_sweep.o
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(MODS)

allclean: clean
	rm -f $(EXE)

.PHONY: all clean allclean
