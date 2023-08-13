rm runlog.dat
touch runlog.dat

cp inputsys-85-130.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-131.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-132.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-133.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-134.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-135.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-136.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-137.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-138.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-139.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-140.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-141.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-142.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-143.dat inputsys.dat
mpiexec -n 11 ./opsplit10v7

cp inputsys-85-144.dat inputsys.dat
exec mpiexec -n 11 ./opsplit10v7

cp runlog.dat runlog-85.dat

