#All source files are places in the SRCS definition with the same
#file names with a .o extension in the OBJS definition

SRCS=oprand85275.f95

OBJS=affall.o tuball.o dvode10.o  tub1d.o tub2d.o tub3d.o tub4d.o tub5d.o tub6d.o tub7d.o tub8d.o tub9d.o tub10d.o oprand85275.o

my_program=oprand85275

CC=mpifort

.SUFFIXES: .o.f95

%.o : %.f95 ; $(CC) -c $*.f95 $(COPTS)

all: $(my_program)

print: 
	lpr $(SRCS)

$(my_program): $(OBJS)
	$(CC) -o $(my_program) $(COPTS) $(OBJS) $(LIBS)