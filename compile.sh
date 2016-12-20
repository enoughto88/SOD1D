rm -r output
mkdir output
cd output
mkdir datafile distribution
cd ..
rm nohup.out
rm -r ifortfile
mkdir ifortfile
ifort -O3 -ipo -parallel -traceback -CB -gen_interfaces -warn all -o SOD.exe \
program/parameter.f90 \
program/euler1d.f90 \
program/main.f90
mv *.mod ifortfile/
mv *__genmod.f90 ifortfile/

