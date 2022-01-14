#!/bin/bash -u

##########################################################################
# Genetic Algorithm                                                      #
# RunCarpSim_V4.sh                                                       #
# This shell script runs openCarp simulations and checks for quiescence  #
# Authors: Kelly Zhang, Chelsea E. Gibbs, Patrick M. Boyle               #
# Version: 4.0.0                                                         #
# Last updated: 1/10/2021                                                #
##########################################################################

# Get commandline parameters
generation=$1
run_type=$2
model_name=$3
experiment=$4
pop_size=$5
sim_dur=$6

# Generate path to new generation directory
NewGenPath="$(pwd)/${experiment}/generation_${generation}";

# Generate path to .par file
if [[ ${run_type} == "test_runs" ]]
then
    # If on test runs, the .par file is unique  within each generation directory
    ParFilePath=${NewGenPath}
else
    # If on full runs, the .par file is overwritten in the experiment directory
    # (All parameters are consistent except for gen  number)
    # Making copies of the full run .par might be extraneous
    ParFilePath="$(pwd)/${experiment}";
fi


# Within the new generation dir, make a dir for run-type (test or full)
# mkdir -p ${NewGenPath}/${run_type};

echo "Checkpoint passed: running carp";
mpiexec -n 6 /Software/cme/openCARP/bin/openCARP.opt +F ${ParFilePath}/${run_type}.par -simID ${NewGenPath}/${run_type} -tend ${sim_dur} -external_imp[0] $(pwd)/${model_name}.so &>${NewGenPath}/${run_type}.log;

echo "Checkpoint passed: running igbextract";
/Software/cme/openCARP/tools/igbutils/igbextract ${NewGenPath}/${run_type}/vm.igb --output-file ${NewGenPath}/${run_type}/vm.nodal --format=asciiTm;

echo "Checkpoint passed: running ExplodeTracesAsciiTm.pl";
perl ExplodeTracesAsciiTm.pl ${NewGenPath}/${run_type}/vm.nodal ${NewGenPath}/${run_type}/vm_exploded;

# TO_KEEP and TO_DISCARD contain the voltage over time data
# for quiescent and deranged variants respectively
TO_KEEP="";
TO_DISCARD="";

# Test_run check for quiescence
if [[ "${run_type}" == "test_runs" ]];
then
    # thresh_time is half of the stim_bcl, which also equals the sim_dur for test runs
    thresh_time=`expr $sim_dur / 2`
    
    # Iterate through the variants in the test pool
    variant_counter=1;
    for FOUT in ${NewGenPath}/${run_type}/vm_exploded*
    do
	
	# Get the most positive voltage after half the sim_dur time and save as variable end_val
        end_val="$(cat ${FOUT} | awk -v thresh_val="$thresh_time" '{if(flag==0||($1>=thresh_val&&$2>max)){flag=1;max=$2;}}END{print(sprintf("%.1f",max));}')"
	
	# If end_val is a more negative voltage than -60 mV then the variant passes the quiescence check
	if [[ "$(echo "${end_val} < -60" | bc)" -eq 1 ]];
	then
	    # echo the variant number into a text file that the .m will read from
	    echo -e "${variant_counter}\n" >>$NewGenPath/${run_type}/quiescent_variants.txt
	    SPONT="${end_val}"
	    TO_KEEP="$TO_KEEP $FOUT";	   
	else
	    SPONT="yes";
	    TO_DISCARD="$TO_DISCARD $FOUT";
	fi
	# Save a log of which variants are spontaneous and quiescent
	echo -e "\n $FOUT $variant_counter $end_val $SPONT" >>./${experiment}/generation_${generation}/spontCheck_log.txt
	variant_counter=$((variant_counter+1))
    done
    
    # Plot the quiescent and deranged traces as sanity check
    xmgrace $TO_KEEP;
    xmgrace $TO_DISCARD;
fi
exit 0;
