#All source files are places in the SRCS definition with the same
#file names with a .o extension in the OBJS definition

SRCS=opsplit10.f95

OBJS=opsplit10.o dvode10.o affall.o tub1.o tub2.o tub3.o tub4.o tub5.o tub6.o tub7.o tub8.o tub9.o tub10.o

my_program=opsplit10

#COPTS=-fp-model precise
CC=mpifort

.SUFFIXES: .o.f95

%.o : %.f95 ; $(CC) -c $*.f95 $(COPTS)

all: $(my_program)

print: 
	lpr $(SRCS)

$(my_program): $(OBJS)
	$(CC) -o $(my_program) $(COPTS) $(OBJS) $(LIBS)
