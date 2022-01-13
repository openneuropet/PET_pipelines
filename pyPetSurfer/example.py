#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Perform regional kinetic modeling (MRTM2) in an example subject using PetSurfer

Note:
    - Path to files and directory have to be modified according to data location
    - This code does not perform the segmentation of the reference region
"""

import os

import numpy as np
import nibabel as nib
import pet_surfer as ps

from importlib import reload
from os.path import join, isfile, isdir
from matplotlib import pyplot as plt

# Freesurfer wrappers from Nipype
from nipype.interfaces.freesurfer.registration import MRICoreg
from nipype.interfaces.freesurfer.model import Binarize
from nipype.interfaces.freesurfer.utils import SampleToSurface
from nipype.interfaces.freesurfer.preprocess import ApplyVolTransform

reload(ps)  # account for any chages to pet_surfer

#%% Download the example dataset from openneuro, if not already available
# https://openneuro.org/datasets/ds001421

# Assume that openneuro is installed and correctly logged in
bids_dir = 'ds001421-download'

if not isdir(bids_dir):
    print('\nDirectory ' + bids_dir + ' does not exist. Attempting to dowload it using openneuro-cli (may take a few minutes)')
    ps.run('openneuro download --snapshot 1.2.1 ds001421 ' + bids_dir)

if not isdir(bids_dir):
    raise ValueError('The example dataset is not available. Please download it manually from https://openneuro.org/datasets/ds001421')

#%% Path to files and directory have to be modified according to data location

# Full path to BIDS PET data
subject = 'sub-01'
session = 'ses-baseline'
base = '_'.join([subject, session])
base_dir = join(bids_dir, subject, session)
pet_file = join(base_dir, 'pet', base + '_pet.nii.gz')
json_file = join(base_dir, 'pet', base + '_pet.json')

# Full path to directory containing FreeSurfer recon (SUBJECTS_DIR)
recon_dir = join(bids_dir, 'derivatives', 'freesurfer')  # i.e., full path to SUBJECTS_DIR

# Create output directory
out_dir = join(bids_dir, 'derivatives', 'pet_surfer')
ps.assert_dir(out_dir)

#%% Create (weighted) average PET images

avg_pet_file = join(out_dir, base + '_pet_weighted_average.nii.gz')

if not isfile(avg_pet_file):
    print('\n--- Creating average PET image ---\n')
    ps.create_weighted_average_pet(
            pet_file,
            json_file,
            avg_pet_file,
            frames=range(9, 26)   # frames 10-26 with 0-based indexing
        )
   
#%% Align PET & MR

lta_file = join(out_dir, base + '_pet_to_anat.lta')
reference_file = join(recon_dir, base, 'mri', 'norm.mgz')

# Coregister average PET & T1
if not isfile(lta_file):
    print('\n--- Aligning PET & MR ---\n')
    coreg = MRICoreg(
            source_file = avg_pet_file,
            reference_file = reference_file,
            out_lta_file = lta_file
        )
    print(coreg.cmdline)
    coreg.run()

# Transform average PET to anatomical space (i.e., Freesurfer recon)

transformed_file = join(out_dir, base + '_space-anat_pet_weighted_average.nii.gz')

if not isfile(transformed_file):
    print('\n--- Transforming average PET to MR space for QC ---\n')
    applyreg = ApplyVolTransform(
            source_file = avg_pet_file,
            target_file = reference_file,
            lta_file = lta_file,
            transformed_file = transformed_file,
            interp = 'cubic'
        )
    print(applyreg.cmdline)
    applyreg.run()

# Visual QC for PET & MR coregistration
qc_png = join(out_dir, 'QC_pet_anat_coreg.png')
if not isfile(qc_png):
    ps.visual_coreg_QC(
            reference_file,
            transformed_file,
            join(out_dir, 'QC_pet_anat_coreg.png'),
            color_range=[50, 800],
            opacity=0.3
        )

#%% Run gtmseg

gtmseg_file = join(recon_dir, base, 'mri', 'gtmseg.mgz')

if not isfile(gtmseg_file):
    
    print('\n--- Computing GTM segmentation ---\n')
    
    gtmseg = ps.GTMSeg(
            subjects_dir = recon_dir,
            subject_id = base,
            out_file = gtmseg_file,
            xcerseg = True,
            no_vermis = True  # do not segment vermis, this is done externally
        )
    print(gtmseg.cmdline)
    gtmseg.run()

# Visual QC for gtmseg
qc_png = join(out_dir, 'QC_gtmseg.png')

if not isfile(qc_png):
    ps.visual_gtmseg_QC(
            join(recon_dir, base),
            join(out_dir, 'QC_gtmseg.png')
        )


#%% Extract volume TACs

# Path to segmentation files - to be modified according to the actual location of the files
cereb_gm_file = join(bids_dir, 'derivatives', 'suit', base + '_space-gtmseg_label-cereb-gm_mask.nii.gz')  # reference region

# The labels_dct dictionary specifies the labels for extracting average TACs.
# The dictionary key is used as output file name. Additional information is 
# specified through the following fields:
#   file: path of the file containing the segmentation
#   ids: labels for which to extract the TACs. List containing int or lists of int.
#        Int will be considered as individual regions, and lists of int specify
#        regions formed by multiple labels. For example,
#        ids: [1, 2, [2, 3], [3, 4, 5]] will output the average TACs for labels
#        1, 2 and the combined region formed by the labels [2, 3] and [3, 4 ,5]. 
#   ext: extension of the output file. Depending on their usage, 
#        "nii.gz" or "dat" have to be used.

labels_dct = {
        base + '_tac-gtmseg-subcort': {
                'file': gtmseg_file,
                # Talamus, caudate, putamen, pallidum, hippocampus, amygdala
                'ids': [10, 11, 12, 13, 17, 18, 49, 50, 51, 52, 53, 54],
                'ext': 'nii.gz'
            },
        base + '_tac-reference': {
                'file': cereb_gm_file,
                'ids': [1],
                'ext': 'dat'
            },
        base + '_tac-high-binding': {
                'file':  gtmseg_file,
                'ids': [[11, 12, 50, 51]],  # indices in the inner list will be considered 
                'ext': 'nii.gz'
            }
    }

# Transform PET data to target space (gtmseg)
reg_pet_file = join(out_dir, base + '_space-gtmseg_pet.nii.gz')

if not ps.tac_file_exist(out_dir, labels_dct):
    
    print('\n--- Extracting Time Activity Curves from volume ---\n')
    
    if not isfile(reg_pet_file):
        print('Transforming dynamic PET file to gtmseg space... (may take a while)')
        applyreg = ApplyVolTransform(
                source_file = pet_file,
                target_file = gtmseg_file,
                lta_file = lta_file,
                transformed_file = reg_pet_file,
                interp = 'nearest'  # trilinear and cubic are otherwise very long to compute
            )
        print(applyreg.cmdline)
        applyreg.run()
    
    # Extract TACs from volume
    ps.extract_vol_tacs(reg_pet_file, out_dir, labels_dct)

#%% Extract regional surface tacs from aparc

for hemi in ['lh', 'rh']:

    fout = join(out_dir, '_'.join([base, 'hemi-' + hemi, 'tac-annot.nii.gz']))
        
    if not isfile(fout):

        print('\n--- Extracting Time Activity Curves from hemi-%s ---\n' % hemi)          

        # Sample TACs to subject surface, no smoothing        
        surf_tacs = join(out_dir, '_'.join([base, 'hemi-' + hemi, 'pet.nii.gz']))
        if not isfile(surf_tacs):
            sampler = SampleToSurface(
                    subjects_dir = recon_dir,
                    subject_id = base,
                    override_reg_subj = True,
                    hemi = hemi,
                    source_file = pet_file,
                    out_file = surf_tacs,
                    reg_file = lta_file,      
                    sampling_method = "point",                
                    sampling_range = 0.5,  # mid pial-wm        
                    sampling_units = "frac",
                    cortex_mask = True
                )
            print(sampler.cmdline)
            sampler.run()
    
        # Extract regional TACs
        annot = join(recon_dir, base, 'label', hemi + '.aparc.annot')
        ps.extract_surf_tacs(surf_tacs, hemi, annot, fout)

#%% Compute MRTM model in high binding regions to obtain an estimate of k2'

ref_file = join(out_dir, base + '_tac-reference.dat')  # TAC of reference region
mrtm_hb_dir = join(out_dir, base + '_high-binding_km-mrtm')  # Output directory
k2p_file = join(mrtm_hb_dir, 'k2prime.dat')

# Create midframe timing file
mid_frames_dat = join(out_dir, 'midframes.sec.dat')
if not isfile(mid_frames_dat):
    ps.create_mid_frame_dat(json_file, mid_frames_dat)

if not isfile(k2p_file):

    print('\n--- Estimating k2prime from high binding regions using MRTM1 ---\n')
    
    # Performing MRTM modeling of high binding region    
    mrtm = ps.MRTM(
            in_file = join(out_dir, base + '_tac-high-binding.nii.gz'),  # TAC of the high-binding region
            mrtm = (ref_file, mid_frames_dat),
            glm_dir = mrtm_hb_dir,
            no_est_fwhm =True
        )
    print(mrtm.cmdline)
    mrtm.run()
    
    # Plot reference and high binding TACs for QC
    timing = np.loadtxt(join(out_dir, 'midframes.sec.dat'), dtype=float)
    ref = np.loadtxt(join(out_dir, base + '_tac-reference.dat'), dtype=float)
    hb = nib.load(join(out_dir, base + '_tac-high-binding.nii.gz')).get_fdata().reshape(-1)
    k2p = np.loadtxt(k2p_file, dtype=float)
    
    plt.plot(timing, np.vstack((ref, hb)).T)
    plt.xlabel('Time (sec)')
    plt.ylabel('Counts')
    plt.title('k2prime: %f' % k2p)
    plt.savefig(join(out_dir, 'ref_hb_tac.png'))

#%% Compute MRTM2 for target regions

# Assign parameterse
tac_files= ['_'.join([base, 'hemi-lh', 'tac-annot']),
            '_'.join([base, 'hemi-rh', 'tac-annot']),
            base + '_tac-gtmseg-subcort']
k2p = np.loadtxt(join(mrtm_hb_dir, 'k2prime.dat'), dtype=float)

# Instantiate nipype MRTM2 model
mrtm2 = ps.MRTM2(
        mrtm2 = (ref_file, mid_frames_dat, k2p),
        no_est_fwhm =True,
        nii_gz=True
    )

# Perform MRTM2 kinetic modeling
for tac_file in tac_files:
    
    glm_dir = glm_dir = join(out_dir, '_'.join([base, tac_file, 'km-mrtm2']))
    bp_file = join(glm_dir, 'bp.nii.gz')
    
    if not isfile(bp_file): 
        print('\n--- Computing regional MRTM2 %s ---\n' % tac_file)
        mrtm2.inputs.in_file = join(out_dir, tac_file + '.nii.gz')  
        mrtm2.inputs.glm_dir = glm_dir  # output directory
        print(mrtm2.cmdline)
        mrtm2.run()

#%% Sample PET data to fsaverage surface, smooth (10 FWHM) and perform 
# vertice-wise kinetic modeling

# Note: dynamic link to the fsaverage directory must be present in recon_dir
# e.g., cd recon_dir; ln -s $FREESURFER_HOME/subject/fsaverage recon_dir/fsaverage

# Make sure fsaverage is in recon_dir
if not isdir(join(recon_dir, 'fsaverage')):
    print('Could not locate fsaverage in ' + recon_dir)
    print('Attempting to create symlink...')
    ps.run('ln -s $SUBJECTS_DIR/fsaverage ' + join(recon_dir, 'fsaverage'))

surf_fwhm = 10
k2p = np.loadtxt(k2p_file, dtype=float)

mrtm2 = ps.MRTM2(
        mrtm2 = (ref_file, mid_frames_dat, k2p),
        no_est_fwhm =True,
        subject_id = 'fsaverage',
        surf = True,
        nii_gz=True
    )

for hemi in ['lh', 'rh']:
    
    glm_dir = join(out_dir, '_'.join([base, 'hemi-' + hemi, 'sm%i' % surf_fwhm, 'km-mrtm2']))
    bp_file = join(glm_dir, 'bp.nii.gz')

    if not isfile(bp_file):
    
        print('\n--- Computing surface (vertice-wise) MRTM2 for hemi-%s ---\n' % hemi)    
    
        # Sample TACs to fsaverage surface and smooth
        
        surf_tacs = join(out_dir, '_'.join([base, 'hemi-' + hemi, 'sm%i' % surf_fwhm, 'pet.nii.gz']))
        
        if not isfile(surf_tacs):
            sampler = SampleToSurface(
                    subjects_dir = recon_dir,
                    subject_id = base,
                    override_reg_subj = True,
                    target_subject = 'fsaverage',
                    hemi = hemi,
                    source_file = pet_file,
                    out_file = surf_tacs,
                    reg_file = lta_file,            
                    sampling_method = "point",             
                    sampling_range = 0.5,  # mid pial-wm        
                    sampling_units = "frac",
                    cortex_mask = True,
                    smooth_surf = surf_fwhm
                )
            print(sampler.cmdline)
            sampler.run()
        
        # Run kinetic modeling (MRTM2)
        mrtm2.inputs.in_file = surf_tacs 
        mrtm2.inputs.glm_dir = glm_dir  # output directory
        mrtm2.inputs.hemi = hemi
        print(mrtm2.cmdline)
        mrtm2.run()
        
#%% Perform voxel-wise kinetic modeling of PET data in MNI152 space

vol_fwhm = 5
glm_dir = join(out_dir, '_'.join([base, 'space-cvs-MNI152', 'sm%i' % vol_fwhm, 'km-mrtm2']))
bp_file = join(glm_dir, 'bp.nii.gz')

if not isfile(bp_file):
    
    print('\n--- Computing volume (voxel-wise) MRTM2 ---\n')
    
    # Transfer PET data to MNI152 space
    # In this specific case, we assume that mri_cvs_register was previously ran:
    # e.g., mri_cvs_register --mov subject --template cvs_avg35_inMNI152
    
    pet_mni152_file = join(out_dir, base + '_space-cvs-MNI152_pet.nii.gz')
    
    if not isfile(pet_mni152_file):
        gcam = ps.GCAM(
                gcam = (
                        pet_file, 
                        lta_file,
                        join(recon_dir, base, 'cvs', 'final_CVSmorph_tocvs_avg35_inMNI152.m3z'),
                        join('$FREESURFER_HOME', 'subjects', 'cvs_avg35_inMNI152', 'mri.2mm', 'register.lta'),
                        '0',
                        1,  # nearest neighbor interpolation
                        pet_mni152_file
                    )
                )
        print(gcam.cmdline)
        gcam.run()
    
    # Create mask
    mask_file = join(out_dir, base + '_space-cvs-MNI152_brain_mask.nii.gz')
    if not isfile(mask_file):
        binvol = Binarize(
                in_file=join(os.environ['FREESURFER_HOME'], 'subjects', 'cvs_avg35_inMNI152', 'mri.2mm', 'aseg.mgz'),
                min=1,
                binary_file=mask_file
            )
        binvol.run()
    
    #% Smooth data
    smooth_pet_file = join(out_dir, base + '_space-cvs-MNI152_sm%i_pet.nii.gz' % vol_fwhm)
    
    if not isfile(smooth_pet_file):
        smoothvol = ps.SmoothVol(
                in_file = pet_mni152_file,
                out_file = smooth_pet_file,
                fwhm = 5,
                mask_file = mask_file
            )
        print(smoothvol.cmdline)
        smoothvol.run()
    
    # Run kinetic modeling (MRTM2)
    mrtm2 = ps.MRTM2(
            mrtm2 = (ref_file, mid_frames_dat, k2p),
            no_est_fwhm =True,
            in_file = smooth_pet_file,
            glm_dir = glm_dir,  # output directory
            mask_file = mask_file,
            nii_gz=True
        )
    print(mrtm2.cmdline)
    mrtm2.run()
    
#%%

print('\nAll computations have been successfully executed.\n')