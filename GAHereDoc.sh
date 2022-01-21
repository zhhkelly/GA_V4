#!/bin/bash -u
run_type=$1
experiment=$2
model_name=$3
generation=$4
pop_size=$5
stim_bcl=$6
stim_num=$7
stim_dur=$8
stim_curr=$9
sim_dur=${10}
test_run_factor=${11}
num_adj=`expr $# - 11`;
electrode_test=`expr ${test_run_factor} \* ${pop_size} \* 1000`
electrode_full=`expr ${pop_size} \* 1000`

if [[ ${run_type} == "test_runs" ]]
then
    ParFile_path="$(pwd)/${experiment}/generation_${generation}"
else
    ParFile_path="$(pwd)/${experiment}"
fi
echo "The ParFile_path is $ParFile_path"

# portion that is the same for all simulations
cat<<START>${ParFile_path}/${run_type}.par
meshname                  =  ${ParFile_path}/${run_type}
simID                     =  ${experiment}

num_phys_regions          =  1
phys_region[0].ptype      =  0
phys_region[0].num_IDs    =  1
phys_region[0].ID         =  0

num_external_imp          =  1
external_imp[0]           = $(pwd)/${model_name}.so

num_imp_regions           =  1

imp_region[0].im          =  ${model_name}
imp_region[0].num_IDs     =  1
imp_region[0].ID          =  0

#Adjustment Files
num_adjustments          =  ${num_adj}

START

counter=0;
while [ "$#" -gt 11 ]; do
cat<<GA_Adjustments>>${ParFile_path}/${run_type}.par
adjustment[${counter}].variable   =  ${model_name}.${12}
adjustment[${counter}].file       =  $(pwd)/${experiment}/generation_${generation}/${run_type}/generation_${generation}_${12}.adj
GA_Adjustments
shift
let "counter+=1"
done;

cat<<MID>>${ParFile_path}/${run_type}.par

##########################################################################################
# Conductivity setup #####################################################################
##########################################################################################
num_gregions              = 1
gregion[0].num_IDs        = 1
gregion[0].ID             = 0
gregion[0].g_mult         = 0

MID

if [[ "${ParFile_path}/${run_type}.par" == "$(pwd)/${experiment}/generation_${generation}/test_runs.par" ]];
then
cat<<TEST>>${ParFile_path}/${run_type}.par
##########################################################################################
# Stimulus settings ######################################################################
##########################################################################################
num_stim		   = 1
stim[0].crct.type          = 0
stim[0].elec.geom_type     = 2
stim[0].elec.p0            = -100 -100 -100
stim[0].elec.p1            = ${electrode_test} 200 200
stim[0].ptcl.start         = 0
stim[0].ptcl.npls          = 0
stim[0].ptcl.duration      = ${stim_dur}
stim[0].pulse.strength     = ${stim_curr}

TEST

elif [[ "${ParFile_path}/${run_type}.par" == "$(pwd)/${experiment}/full_runs.par" ]];
then
cat<<FULL>>${ParFile_path}/${run_type}.par
##########################################################################################
# Stimulus settings ######################################################################
##########################################################################################
num_stim		   = 1
stim[0].crct.type          = 0
stim[0].elec.geom_type     = 2
stim[0].elec.p0            = -100 -100 -100
stim[0].elec.p1            = ${electrode_full} 200 200
stim[0].ptcl.start         = 0
stim[0].ptcl.duration      = ${stim_dur}
stim[0].pulse.strength     = ${stim_curr}
FULL
if ( [ ${stim_num} == 0 ] && [ ${stim_bcl}==0 ] );
then
cat<<NOSTIM>>${ParFile_path}/${run_type}.par

NOSTIM
else
cat<<STIM>>${ParFile_path}/${run_type}.par
stim[0].ptcl.npls          = ${stim_num}
stim[0].ptcl.bcl           = ${stim_bcl}

STIM
fi
fi
cat<<EOF>>${ParFile_path}/${run_type}.par
##########################################################################################
# openCARP settings ######################################################################
##########################################################################################
tend                      =  ${sim_dur}
timedt                    =  10
spacedt                   =  .1
dt                        =  10

num_LATs                  =  1
lats[0].measurand         =  0
lats[0].all               =  1
lats[0].threshold         =  0.0
lats[0].method            =  1

experiment                =  0
bidomain                  =  0
mass_lumping              =  0
parab_solve               =  1

EOF

echo "Checkpoint passed: Made .par file for ${run_type}"
