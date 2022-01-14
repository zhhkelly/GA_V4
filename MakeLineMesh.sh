#!/bin/bash -u
##########################################################################
# Genetic Algorithm                                                      #
# MakeLineMesh.sh                                                        #
# This shell script creates the line mesh files for running a population #
# Authors: Kelly Zhang, Chelsea E. Gibbs, Patrick M. Boyle               #
# Version: 4.0.0                                                         #
# Last updated: 1/10/2021                                                #
##########################################################################

mesh_path=$1
N=$2;
MESH=$3;


echo $N > $mesh_path/$MESH.pts;
for N in $(seq 1 $N); do
    echo $N |awk '{print(100*$1,0,0);}' >> $mesh_path/$MESH.pts;
done;

echo $((N-1)) > $mesh_path/$MESH.elem;
echo "1" > $mesh_path/$MESH.lon;
for N in $(seq 2 $N); do
    echo $N |awk '{print("Ln", $1-2, $1-1, 0);}' >> $mesh_path/$MESH.elem;
    echo "1 0 0" >> $mesh_path/$MESH.lon;
done;
echo "Checkpoint passed: ${MESH} mesh made for ${mesh_path}"
exit 0;
