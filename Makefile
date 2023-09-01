#All source files are places in the SRCS definition with the same
#file names with a .o extension in the OBJS definition

SRCS=opsplit10v7.f90

OBJS=opsplit10v7.o  affall.o tuball.o dvode10.o tub1.o tub2.o tub3.o tub4.o tub5.o tub6.o tub7.o tub8.o tub9.o tub10.o

my_program=opsplit10v7

#COPTS=-fp-model precise
CC=mpifort

.SUFFIXES: .o.f90

%.o : %.f90 ; $(CC) -c $*.f90 $(COPTS)

all: $(my_program)

print: 
	lpr $(SRCS)

$(my_program): $(OBJS)
	$(CC) -o $(my_program) $(COPTS) $(OBJS) $(LIBS)
