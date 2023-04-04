CC       = gcc
CCOPTS   = -c -O
LINK     = gcc
LINKOPTS = liblpsolve55.a -lm -ldl

.c.o: 
	$(CC) $(CCOPTS) $<

all:  computeflowsingle  

computeflow.o: computeflow.c

computeflowsingle.o: computeflowsingle.c


computeflowsingle: computeflow.o computeflowsingle.o
	$(LINK) -o computeflowsingle computeflowsingle.o computeflow.o $(LINKOPTS)

clean:
	rm *o computeflowsingle  

