# template for the Intel fortran compiler
# typical use with mkmf
# mkmf -t template.ifc -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include
FC = ifc
LD = ifc
FFLAGS = -i4 -r8 -fpp -O2 -I/opt/Intel/include
# FFLAGS = -i4 -r8 -fpp -g #debugging
LDFLAGS = $(LIBS)
# LIBS needs to be customized for your site
LIBS = -L/usr/local/lib -lPEPCF90 -L/opt/Intel/lib -ludunits -lnetcdf -L/usr/mpi/lib
