# debug
gfortran -c -fdefault-real-8 -fbackslash -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace ray_bilinear_face_intersection.f90 
gfortran -fdefault-real-8 -fbackslash -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace ray_bilinear_face_intersection.o -o raybi
## optimized
#gfortran -c -O3 ray_bilinear_face_intersection.f90 
#gfortran -O3 ray_bilinear_face_intersection.o -o raybi
## cleaning 

# debug
#gfortran -c -fdefault-real-8 -fbackslash -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace ray_speedtest.f90 
#gfortran -fdefault-real-8 -fbackslash -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace    ray_speedtest.o -o raybi
# optimized
#gfortran -c -fdefault-real-8 -O3 ray_speedtest.f90 
#gfortran    -fdefault-real-8 -O3 ray_speedtest.o -o raybi
# cleaning 

rm *.o 
echo ' Compilation done.'
