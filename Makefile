.SUFFIXES :
.SUFFIXES : .f90 .o

FC = gfortran
#FFLAGS = -O2 -I/opt/local/include -fbounds-check -ftrace=full
FFLAGS = -O2 -I/usr/local/include -Wall -g -fcheck=array-temps,bounds,do,mem,pointer,recursion -ffpe-trap=invalid,zero,overflow,underflow
LD = f90
LDFLAGS = -L/usr/local/lib
LDLIBS = -lfftw3

SRC = math_module.f90 planet_module.f90 \
	fft_module.f90 glatwgt_module.f90 alf_module.f90 \
	grid_module.f90 time_module.f90 uv_module.f90 \
	legendre_transform_module.f90 init_module.f90 \
	upstream_module.f90 interpolate_module.f90 \
	polint_module.f90 bicubic_module.f90 \
	euler_module.f90 semilag_module.f90 nisl_module.f90 \
	sphere_module.f90  field_module.f90 main.f90
OBJ = ${SRC:.f90=.o}
TARGET=adv

$(TARGET) : $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) $(LDLIBS) -o $@

time_module.o: grid_module.o planet_module.o
glatwgt_module.o : math_module.o
legendre_transform_module.o : glatwgt_module.o alf_module.o fft_module.o 
init_module.o : math_module.o planet_module.o sphere_module.o legendre_transform_module.o uv_module.o
upstream_module.o : sphere_module.o grid_module.o time_module.o interpolate_module.o
interpolate_module.o : math_module.o grid_module.o sphere_module.o bicubic_module.o polint_module.o
grid_module.o : legendre_transform_module.o init_module.o uv_module.o
euler_module.o : planet_module.o grid_module.o time_module.o legendre_transform_module.o uv_module.o
semilag_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_module.o 
nisl_module.o : grid_module.o time_module.o legendre_transform_module.o upstream_module.o sphere_module.o field_module.o
main.o : grid_module.o time_module.o euler_module.o semilag_module.o nisl_module.o field_module.o

clean :
	rm -f *.o *.mod $(TARGET) *.dat $(TARGET).log

.f90.o :
	$(FC) $(FFLAGS) $< -c
