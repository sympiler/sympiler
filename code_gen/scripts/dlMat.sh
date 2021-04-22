#!/bin/sh
OUTPUT=$1
binFile=$2
#This script downlds the SC17 matrix set and extract it to output folder and then generates .
wget http://www.cise.ufl.edu/research/sparse/MM/Oberwolfach/gyro_k.tar.gz -O ${OUTPUT}gyro_k.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/UTEP/Dubcova2.tar.gz -O ${OUTPUT}Dubcova2.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Boeing/msc23052.tar.gz -O ${OUTPUT}msc23052.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/ACUSIM/Pres_Poisson.tar.gz -O ${OUTPUT}Pres_Poisson.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/TKK/cbuckle.tar.gz -O ${OUTPUT}cbuckle.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Botonakis/thermomech_dM.tar.gz -O ${OUTPUT}thermomech_dM.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Simon/olafu.tar.gz -O ${OUTPUT}olafu.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/UTEP/Dubcova3.tar.gz -O ${OUTPUT}Dubcova3.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Wissgott/parabolic_fem.tar.gz -O ${OUTPUT}parabolic_fem.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/McRae/ecology2.tar.gz -O ${OUTPUT}ecology2.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Oberwolfach/gyro.tar.gz -O ${OUTPUT}gyro.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Simon/raefsky4.tar.gz -O ${OUTPUT}raefsky4.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/CEMW/tmt_sym.tar.gz -O ${OUTPUT}tmt_sym.tar.gz
wget http://www.cise.ufl.edu/research/sparse/MM/Rajat/rajat21.tar.gz -O ${OUTPUT}rajat21.tar.gz


for f in ${OUTPUT}*.gz; do tar xvf $f -C ${OUTPUT}; done;
#copying all matrices into root so it is easier to explore later.
cp ${OUTPUT}gyro_k/gyro_k.mtx ${OUTPUT}gyro_k.mtx;
cp ${OUTPUT}Dubcova2/Dubcova2.mtx ${OUTPUT}Dubcova2.mtx;
cp ${OUTPUT}msc23052/msc23052.mtx ${OUTPUT}msc23052.mtx;
cp ${OUTPUT}Pres_Poisson/Pres_Poisson.mtx ${OUTPUT}Pres_Poisson.mtx;
cp ${OUTPUT}cbuckle/cbuckle.mtx ${OUTPUT}cbuckle.mtx;
cp ${OUTPUT}thermomech_dM/thermomech_dM.mtx ${OUTPUT}thermomech_dM.mtx;
cp ${OUTPUT}olafu/olafu.mtx ${OUTPUT}olafu.mtx;
cp ${OUTPUT}Dubcova3/Dubcova3.mtx ${OUTPUT}Dubcova3.mtx;
cp ${OUTPUT}parabolic_fem/parabolic_fem.mtx ${OUTPUT}parabolic_fem.mtx;
cp ${OUTPUT}ecology2/ecology2.mtx ${OUTPUT}ecology2.mtx;
cp ${OUTPUT}gyro/gyro.mtx ${OUTPUT}gyro.mtx;
cp ${OUTPUT}raefsky4/raefsky4.mtx ${OUTPUT}raefsky4.mtx;
cp ${OUTPUT}tmt_sym/tmt_sym.mtx ${OUTPUT}tmt_sym.mtx;
mkdir ${OUTPUT}ccache
cp ${OUTPUT}rajat21/rajat21.mtx ${OUTPUT}ccache/rajat21.mtx;

#Generating triangular matrices
#mkdir ${OUTPUT}triangular
#for f in ${OUTPUT}*.mtx
#do
#    echo "Generating Triangular matrix for $f"
#    ${binFile} $f ${OUTPUT}triangular/
#done
