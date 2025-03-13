touch runlog.dat

cp inputsys-85-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-90-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-95-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-100-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-105-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-110-321.dat inputsys.dat
mpiexec -n 11 ./arcopmurvar275

cp inputsys-115-321.dat inputsys.dat
exec mpiexec -n 11 ./arcopmurvar275





