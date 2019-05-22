# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -O3

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./ -I../include

# define the library path
LFLAGS	= -L/usr/local/lib -L../lib

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas -lds -lal

# define the C object files
OBJS	= scheme.o

#define the directory for object
OBJSDIR = object

# define the executable file
MAIN	= exe

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS) $(LFLAGS) $(INCLUDES) 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^

lib: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^



# clean the executable file and object files
clean:
	$(RM) $(OBJS) $(MAIN)
