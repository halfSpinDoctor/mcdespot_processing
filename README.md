# mcdespot_processing
Processing script for analysis of mcDESPOT data using the FSL and QUIT software
packages

Pre-processing steps (DICOM->NIfTI conversion, brain masking, registration, and
fitting of T1/B1 maps) are computed on the local computer. 

The computation of T2 and off-resonance maps (DESPOT2-FM) and multi-component
relaxation/fraction maps (mcDESPOT) are submitted as a slurm job using sbatch,
so that they will run on a cluster.

The job is split up by slices (z). This software includes a modified perl
script borrowed from the ANTs software package that will wait for the cluster
job to finish (waitForSlurmJobs_sah.pl), then re-combine the individual slices
into a single volumetric image.

The script takes as an input four files that point to the folders with the DICOM
series for your images:

spgr_list
irspgr_list
ssfp_0_list
ssfp_180_list

Here, you want to provide relative or absolute path to the DICOM folders. A bit
of a gotcha, you have to make sure that there is a return/newline after the last
entry. I.e. the file ends with a blank line at the bottom, otherwise it will
skip reading in the last dataset. Also, as of right now, the flip angles and TR
values are hard-coded, not read from the DICOM header. So if there is some
variation between datasets, or you change the protocol, the script will need to
be updated by hand.

The script will create a number of directories to organise the results of
pre-processing (mask, registration, etc), and to collect the output of the
cluster slurm jobs. These are created from the folder where the script is run:

originalData -
  the original images, in NIfTI format (converted from DICOM using dcm2niix)
  
maskedData -
  brain mask and refweight (used for registration). Images reoriented to axial
  
registredData -
  all images aligned to the 1st SPGR image. Uses flirt with -refweight

singleComponent -
  single component T1 and T2 maps, and B1 and f0 (off-resonance) calibration maps
  
multiComponent -
  output of mcDESPOT multicomponent fitting. See QUIT documentation for details
