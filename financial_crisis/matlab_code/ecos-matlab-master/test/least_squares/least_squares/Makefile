ECOS_PATH = ../../..
INCLUDES = -I$(ECOS_PATH)/code/include -I$(ECOS_PATH)/code/external/SuiteSparse_config
LIBS = -lm $(ECOS_PATH)/code/libecos.a $(ECOS_PATH)/code/external/amd/libamd.a $(ECOS_PATH)/code/external/ldl/libldl.a

all: solver.o
	gcc -Wall -O3 -o testsolver testsolver.c solver.o $(LIBS) $(INCLUDES)

solver.o: solver.c
	gcc -ansi -Wall -O3 -c solver.c $(INCLUDES)

clean:
	rm testsolver *.o
