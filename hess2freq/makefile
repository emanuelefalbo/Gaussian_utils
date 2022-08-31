#Start of Makefile
FC=gfortran
FFLAGS= -O3 -Wall -Wextra 
LIB= -llapack -lblas
SRC= read_hess.f95 \
     hess_main.f95   
OBJ=${SRC:.f95=.o}

%.o: %.f95
	 $(FC) $(FFLAGS) -o $@ -c $< 
	 
hess2freq.exe: $(OBJ)
		 $(FC) $(FFLAGS) -o $@ $(OBJ) -llapack -lblas

clean:
		@rm -f *.mod *.o hess2freq.exe
