export ALGENCAN=$HOME/algencan-3.1.1

gfortran -c -O3 covering-disconnected.f90	
gfortran -c -O3 covering-heart.f90		
gfortran -c -O3 covering-rings.f90		
gfortran -c -O3 covering-soap.f90 	
gfortran -c -O3 covering-star.f90 	
gfortran -c -O3 covering-twosq.f90	
gfortran -c -O3 covering.f90              

gfortran covering-disconnected.o -L$ALGENCAN/lib -lalgencan -o covering-disconnected 
gfortran covering-heart.o        -L$ALGENCAN/lib -lalgencan -o covering-heart	     
gfortran covering-rings.o        -L$ALGENCAN/lib -lalgencan -o covering-rings	     
gfortran covering-soap.o         -L$ALGENCAN/lib -lalgencan -o covering-soap	     
gfortran covering-star.o         -L$ALGENCAN/lib -lalgencan -o covering-star	     
gfortran covering-twosq.o        -L$ALGENCAN/lib -lalgencan -o covering-twosq

./covering-disconnected > covering-disconnected.out
./covering-heart        > covering-heart.out	     
./covering-rings        > covering-rings.out	     
./covering-soap         > covering-soap.out	     
./covering-star         > covering-star.out	     
./covering-twosq        > covering-twosq.out        
