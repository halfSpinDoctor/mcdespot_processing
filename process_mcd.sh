#!/bin/bash

## process_mcd.sh
#
# Shell script to process Karen's mcDESPOT study data using QUIT pacakage.
# mcDESPOT processing is parallelised by slice
#

## BSD 2-Clause License
#
# Copyright (c) 2017, Samuel A. Hurley
# University of Wisconsin - Madison
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## USER SETTINGS
#
# Set scan and fitting parameters here
#
#
# SPGR Scan Parameters
SPGR_FLIP="3 4 5 6 7 9 13 18"
SPGR_TR="0.0053"

# IRSPGR Scan Parameters
IRSPGR_FLIP="5"
IRSPGR_TR="0.0053"
IRSPGR_NPE="96"
IRSPGR_TI="0.450"
#
# SSFP Scan Parameters
SSFP_FLIP="9.7059 12.9412 16.9853 21.8382 26.6912 32.3529 41.2500 55.0000"
SSFP_PHASE="0 180"
SSFP_TR="0.0054"
#
# Set number of threads used for DESPOT2-FM and mcDESPOT (Default is 4)
NUM_THREADS=8
#
# END OF USER SETTINGS


## QUIT SETTINGS
#
# Settings passed to the QUIT software package
#
export QUIT_EXT=NIFTI_GZ
export FSLOUTPUTTYPE=NIFTI_GZ
#
# END OF QUIT SETTINGS

# trap keyboard interrupt (control-c)
trap control_c SIGINT

control_c()
# run if user hits control-c
{
  echo -en "\n*** User pressed CTRL + C ***\n"
  scancel $jobIDs

  exit $?
  echo -en "\n*** Script cancelled by user ***\n"
}

## 1. Create directory structure
mkdir -p originalData/spgr
mkdir -p originalData/irspgr
mkdir -p originalData/ssfp_0
mkdir -p originalData/ssfp_180

mkdir -p maskedData/spgr
mkdir -p maskedData/irspgr
mkdir -p maskedData/ssfp_0
mkdir -p maskedData/ssfp_180

mkdir -p registeredData/spgr
mkdir -p registeredData/irspgr
mkdir -p registeredData/ssfp_0
mkdir -p registeredData/ssfp_180

# output of slurm cluster jobs
mkdir -p slurmOutput

# osingle and multicomponent maps
mkdir singleComponent
mkdir multiComponent

## 2. Convert DICOM-NIFTI

# get number input images from list files
NUM_SPGR=`wc -l < spgr_list`
NUM_IRSPGR=`wc -l < irspgr_list`
NUM_SSFP_0=`wc -l < ssfp_0_list`
NUM_SSFP_180=`wc -l < ssfp_180_list`

# Iterate over input file lists
ii=1;
cat spgr_list | while read fname
do
	echo $fname;
	dcm2niix -b y -z y -f spgr_$ii -o originalData/spgr $fname
	ii=`expr $ii + 1`
done

ii=1;
cat irspgr_list | while read fname
do
	echo $fname;
	dcm2niix -b y -z y -f irspgr_$ii -o originalData/irspgr $fname
	ii=`expr $ii + 1`
done

ii=1;
cat ssfp_0_list | while read fname
do
	echo $fname;
	dcm2niix -b y -z y -f ssfp_0_$ii -o originalData/ssfp_0 $fname
	ii=`expr $ii + 1`
done

ii=1;
cat ssfp_180_list | while read fname
do
	echo $fname;
	dcm2niix -b y -z y -f ssfp_180_$ii -o originalData/ssfp_180 $fname
	ii=`expr $ii + 1`
done

## 3. Create brain mask
#
# FSL BET does not seem to handle the higher flip angle images well, so run bet to generate brain mask from first image,
# The mask will be input to QUIT processing to select which voxels to process
#
# Using a weighting mask with FLIRT registration instead of cropping the image improves the quality of alignment significantly:
# Using the -refweight option does not create artificial boundaries in the image for the registration to pick up on
#

# create brain mask from first spgr
i=1;
	echo Brain extracting spgr_$i
	fslreorient2std originalData/spgr/spgr_$i maskedData/spgr/spgr_std_$i 
	
	# create brain mask 
	bet maskedData/spgr/spgr_std_$i maskedData/spgr/spgr_std_brain_$i -m
	cp maskedData/spgr/spgr_std_brain_"$i"_mask.nii.gz maskedData/brainmask.nii.gz
	
	# create refweight for registration
	fslmaths maskedData/spgr/spgr_std_brain_"$i"_mask.nii.gz -mul 0.5 -add 0.5 maskedData/refweight

## 4. Reorient and brain extract images
#
# Reorient images to axial orientation prior to registration
#

# spgr
for i in `seq 1 1 $NUM_SPGR`
do
	echo Reorienting spgr_$i to axial
	fslreorient2std originalData/spgr/spgr_$i maskedData/spgr/spgr_std_$i 
	# fslmaths maskedData/spgr/spgr_std_$i -mul maskedData/spgr/spgr_std_brain_1_mask maskedData/spgr/spgr_std_brain_$i
done

# irspgr
for i in `seq 1 1 $NUM_IRSGPR`
do
	echo Reorienting irspgr_$i to axial
	fslreorient2std originalData/irspgr/irspgr_$i maskedData/irspgr/irspgr_std_$i
	# fslmaths maskedData/irspgr/irspgr_std_$i -mul maskedData/irspgr/irspgr_std_brain_1_mask maskedData/irspgr/irspgr_std_brain_$i
done

# ssfp_0
for i in `seq 1 1 $NUM_SSFP_0`
do
	echo Reorienting ssfp_0_$i to axial
	fslreorient2std originalData/ssfp_0/ssfp_0_$i maskedData/ssfp_0/ssfp_0_std_$i
	# fslmaths maskedData/ssfp_0/ssfp_0_std_$i -mul maskedData/ssfp_0/ssfp_0_std_brain_1_mask maskedData/ssfp_0/ssfp_0_std_brain_$i
done

# ssfp_180
for i in `seq 1 1 $NUM_SSFP_180`
do
	echo Reorienting ssfp_180_$i to axial
	fslreorient2std originalData/ssfp_180/ssfp_180_$i maskedData/ssfp_180/ssfp_180_std_$i
	# fslmaths maskedData/ssfp_180/ssfp_180_std_$i -mul maskedData/ssfp_180/ssfp_180_std_brain_1_mask maskedData/ssfp_180/ssfp_180_std_brain_$i
done


## 5. Register images
#
#

# for spgr, register images to first flip angle
FLIRTOPTS="-ref maskedData/spgr/spgr_std_1 -dof 6 -nosearch -refweight maskedData/refweight -interp sinc"

# spgr
for i in `seq 1 1 $NUM_SPGR`
do
	echo Registration of spgr_$i
	flirt -in maskedData/spgr/spgr_std_$i -out registeredData/spgr/spgr_$i $FLIRTOPTS
done

# irspgr
for i in `seq 1 1 $NUM_IRSPGR`
do
	echo Registration of irspgr_1
	flirt -in maskedData/irspgr/irspgr_std_$i -out registeredData/irspgr/irspgr_$i $FLIRTOPTS
done

# ssfp_0
for i in `seq 1 1 $NUM_SSFP_0`
do
	echo Registration of ssfp_0_$i
	flirt -in maskedData/ssfp_0/ssfp_0_std_$i -out registeredData/ssfp_0/ssfp_0_$i $FLIRTOPTS
done

# ssfp_180
for i in `seq 1 1 $NUM_SSFP_180`
do
	echo Registration of ssfp_180_$i
	flirt -in maskedData/ssfp_180/ssfp_180_std_$i -out registeredData/ssfp_180/ssfp_180_$i $FLIRTOPTS
done

# Merge outputs
fslmerge -t registeredData/spgr_merge registeredData/spgr/spgr*
fslmerge -t registeredData/irspgr_merge registeredData/irspgr/irspgr*
fslmerge -t registeredData/ssfp_merge registeredData/ssfp_0/ssfp_0* registeredData/ssfp_180/ssfp_180_*

## 6. Compute T1, B1 maps from DESPOT1-HIFI
#
#
# Enter flip-angles (degrees): 3 4 5 6 7 9 13 18
# Enter TR (seconds): 0.0053
# Enter read-out flip-angle (degrees): 5
# Enter read-out TR (seconds): 0.0053
# Enter number of spatial locations (remember +4): 96
# Enter TIs (seconds): 0.450
#

qidespot1hifi --out singleComponent/ --mask maskedData/brainmask.nii.gz --verbose registeredData/spgr_merge.nii.gz registeredData/irspgr_merge.nii.gz << EOL
$SPGR_FLIP
$SPGR_TR
$IRSPGR_FLIP
$IRSPGR_TR
$IRSPGR_NPE
$IRSPGR_TI
EOL


## 7. Compute T2, off-resonance (f0) maps from DESPOT2-FM
#
#
# Enter flip-angles (degrees): 9.7059 12.9412 16.9853 21.8382 26.6912 32.3529 41.2500 55.0000
# Enter phase-increments (degrees): 0 180
# Enter TR (seconds): 0.0054
#

# Generate a bash script to submit a slurm job
# <start: slurm script>
cat > despot2fm_job.sh << EOF
#!/bin/bash

qidespot2fm --start \$SLURM_ARRAY_TASK_ID --stop \`expr \$SLURM_ARRAY_TASK_ID + 1\` --out singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_ --mask maskedData/brainmask.nii.gz --B1 singleComponent/HIFI_B1.nii --threads $NUM_THREADS --verbose singleComponent/HIFI_T1.nii registeredData/ssfp_merge.nii.gz << EOL
$SSFP_FLIP
$SSFP_PHASE
$SSFP_TR
EOL

# Trim the output to a single slice
fslroi singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_f0 singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_f0 \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_PD singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_PD \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_residual singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_residual \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_T2 singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_T2 \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1


# Clean up full volume images (extension of .nii)
rm -f singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_f0 singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_*.nii
EOF
# <end: slurm script>

# Get number of slices for slurm parallelisation [0..(nslice-1)]
nslice=`fslval registeredData/ssfp_merge dim3`
nslice=`expr $nslice - 1`

# Run jobs on SLURM and wait to finish
# Run as 1 task with 8 threads (CPUs) per task. See https://www.massive.org.au/userguide/running-slurm-jobs/running-multi-threading-jobs
jobIDs=`sbatch --array=0-$nslice --ntasks=1 --cpus-per-task=$NUM_THREADS despot2fm_job.sh  | rev | cut -f1 -d\ | rev`
echo
echo "--------------------------------------------------------------------------------------"
echo " Starting DESPOT2-FM on SLURM cluster. Submitted $nslice jobs "
echo "--------------------------------------------------------------------------------------"
    # now wait for the jobs to finish.
	./waitForSlurmJobs_sah.pl 1 60 $jobIDs

# Returns 1 if there are errors
if [ ! $? -eq 0 ]; then
    echo "SLURM submission failed - jobs went into error state"
    mv slurm-*.out slurmOutput
    exit 1;
fi

# Cleanup slurm command line output
mv slurm-*.out slurmOutput

# Merge DESPOT2-FM output slices into single volume
fslmerge -z singleComponent/FM_f0 singleComponent/slice_*_FM_f0.nii.gz
fslmerge -z singleComponent/FM_PD singleComponent/slice_*_FM_PD.nii.gz
fslmerge -z singleComponent/FM_residual singleComponent/slice_*_FM_residual.nii.gz
fslmerge -z singleComponent/FM_T2 singleComponent/slice_*_FM_T2.nii.gz
rm -f singleComponent/slice_*


## 8. Compute multi-component maps from mcDESPOT
#
#
# Enter input filename: registeredData/spgr_merge.nii.gz  
# Reading file: registeredData/spgr_merge.nii.gz
# Enter sequence type (SPGR/SSFP): SPGR
# Enter flip-angles (degrees): 3 4 5 6 7 9 13 18
# Enter TR (seconds): 0.0053
# Enter next filename (END to finish input): registeredData/ssfp_merge.nii.gz
# Reading file: registeredData/ssfp_merge.nii.gz
# Enter sequence type (SPGR/SSFP): SSFP
# Enter flip-angles (degrees): 9.7059 12.9412 16.9853 21.8382 26.6912 32.3529 41.2500 55.0000
# Enter phase-increments (degrees): 0 180
# Enter TR (seconds): 0.0054
# Enter next filename (END to finish input): END

# Normalize input data by proton density maps - NOT USED
fslmaths registeredData/spgr_merge.nii.gz -div singleComponent/HIFI_PD.nii registeredData/spgr_merge_norm
fslmaths registeredData/ssfp_merge.nii.gz -div singleComponent/FM_PD.nii.gz registeredData/ssfp_merge_norm


# Generate a bash script to submit a slurm job
# <start: slurm script>
cat > mcdespot_job.sh << EOF
#!/bin/bash

# qimcdespot --start \$SLURM_ARRAY_TASK_ID --stop \`expr \$SLURM_ARRAY_TASK_ID + 1\` --out multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_ --mask maskedData/brainmask.nii.gz --model 3 --tesla 3 --B1 singleComponent/HIFI_B1.nii --f0 singleComponent/FM_f0.nii.gz --threads $NUM_THREADS --verbose << EOL
# normalise signals to mean (--scale)
qimcdespot --start \$SLURM_ARRAY_TASK_ID --stop \`expr \$SLURM_ARRAY_TASK_ID + 1\` --out multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_ --scale --mask maskedData/brainmask.nii.gz --model 3 --tesla 3 --B1 singleComponent/HIFI_B1.nii --f0 singleComponent/FM_f0.nii.gz --threads $NUM_THREADS --verbose << EOL
registeredData/spgr_merge.nii.gz
SPGR
$SPGR_FLIP
$SPGR_TR
registeredData/ssfp_merge.nii.gz
SSFP
$SSFP_FLIP
$SSFP_PHASE
$SSFP_TR
END
EOL

# Trim the output to a single slice
fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_B1 multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_B1 \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f0 multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f0 \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f_csf multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f_csf \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f_m multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_f_m \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_iterations multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_iterations \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_PD multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_PD \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_csf multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_csf \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_ie multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_ie \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_m multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T1_m \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_csf multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_csf \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_ie multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_ie \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_m multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_T2_m \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1

fslroi multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_tau_m multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_3C_tau_m \
	0 -1 0 -1 \$SLURM_ARRAY_TASK_ID 1


# Clean up full volume images (extension of .nii)
# Not needed if output type is set to NIFTI_GZ
# rm -f singleComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_FM_f0 multiComponent/slice_"\$(printf %04d \$SLURM_ARRAY_TASK_ID)"_*.nii

EOF
# <end: slurm script>

# Get number of slices for slurm parallelisation [0..(nslice-1)]
nslice=`fslval registeredData/ssfp_merge dim3`
nslice=`expr $nslice - 1`

# Run jobs on SLURM and wait to finish
jobIDs=`sbatch --array=0-$nslice mcdespot_job.sh  | rev | cut -f1 -d\ | rev`
# jobIDs=`sbatch --array=0-$nslice --ntasks=1 --cpus-per-task=$NUM_THREADS mcdespot_job.sh  | rev | cut -f1 -d\ | rev`
echo
echo "--------------------------------------------------------------------------------------"
echo " Starting mcDESPOT on SLURM cluster. Submitted $nslice jobs "
echo "--------------------------------------------------------------------------------------"

# now wait for the jobs to finish.
./waitForSlurmJobs_sah.pl 1 60 $jobIDs

# Returns 1 if there are errors
if [ ! $? -eq 0 ]; then
    echo "SLURM submission failed - jobs went into error state"
    mv slurm-*.out slurmOutput
    exit 1;
fi

# Cleanup slurm command line output
mv slurm-*.out slurmOutput

# Merge mcDESPOT output slices into single volume

# Calibration maps (B1, off-resonance [f0])
fslmerge -z multiComponent/3C_B1 multiComponent/slice_*_3C_B1.nii.gz
fslmerge -z multiComponent/3C_f0 multiComponent/slice_*_3C_f0.nii.gz

# Compartment fractions
fslmerge -z multiComponent/3C_f_csf multiComponent/slice_*_3C_f_csf.nii.gz
fslmerge -z multiComponent/3C_f_m multiComponent/slice_*_3C_f_m.nii.gz

fslmerge -z multiComponent/3C_iterations multiComponent/slice_*_3C_iterations.nii.gz
fslmerge -z multiComponent/3C_PD multiComponent/slice_*_3C_PD.nii.gz

# Compartment T1 times
fslmerge -z multiComponent/3C_T1_csf multiComponent/slice_*_3C_T1_csf.nii.gz
fslmerge -z multiComponent/3C_T1_ie multiComponent/slice_*_3C_T1_ie.nii.gz
fslmerge -z multiComponent/3C_T1_m multiComponent/slice_*_3C_T1_m.nii.gz

# Compartment T2 times
fslmerge -z multiComponent/3C_T2_csf multiComponent/slice_*_3C_T2_csf.nii.gz
fslmerge -z multiComponent/3C_T2_ie multiComponent/slice_*_3C_T2_ie.nii.gz
fslmerge -z multiComponent/3C_T2_m multiComponent/slice_*_3C_T2_m.nii.gz

# Myelin resonance time
fslmerge -z multiComponent/3C_tau_m multiComponent/slice_*_3C_tau_m.nii.gz

# Fitting bounds
cat multiComponent/slice_*3C_bounds.txt > multiComponent/3C_bounds.txt

# Clean up
rm -f multiComponent/slice_*
