#! /usr/bin/env python

# -*- coding: utf-8 -*-
"""PetSurfer

Collection of functions and wrappers for performing the analysis of PET data
using FreeSurfer's PET pipeline (https://surfer.nmr.mgh.harvard.edu/fswiki/PetSurfer)
"""

import os
import nibabel as nib
import numpy as np
import json

from subprocess import Popen, PIPE
from os.path import join, isfile

from nipype.interfaces.base import (
    TraitedSpec,
    File,
    traits,
    InputMultiPath,
    OutputMultiPath,
    Directory,
    isdefined,
)
from nipype.interfaces.freesurfer.base import FSCommand, FSTraitedSpec
from nipype.interfaces.freesurfer.preprocess import (
    ApplyVolTransformInputSpec,
    ApplyVolTransformOutputSpec,
    ApplyVolTransform
)

# For MRTM/MRTM2
from nipype.interfaces.freesurfer.model import (
        GLMFitInputSpec,
        GLMFitInputSpec,
        GLMFit
    )

#%% Utility functions


def assert_dir(dir_path):
    
    """
    Create directory, if it does not exist
    
    Arguments
    ---------
    dir_path: string
        path to directory to be create   
    """ 
    
    full_path = os.path.abspath(dir_path)
    if not os.path.isdir(full_path):
        print('Creating %s' % full_path)
        os.makedirs(full_path)


def run(cmd):
    
    """
    Excute a command
    
    Arguments
    ---------
    cmd: string
        command to be executed
    """ 
    
    print('\n' + cmd)
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    output, error = p.communicate()
    if output:
        print(output.decode('latin_1'))
    if error:
        print(error.decode('latin_1'))


#%% Nipype wrappers


class SmoothVolInputSpec(FSTraitedSpec):
    
    in_file = File(
        argstr="--i %s",
        desc="input volume",
        mandatory=True
    )
    
    out_file = File(
        argstr="--o %s",
        desc="save input after smoothing",
    )
    
    save_detrended = traits.Bool(
        argsstr='--save-detrended',
        desc='detrend output when saving'
    )

    save_unmask = traits.Bool(
        argstr='--save-unmasked',
        desc='do not mask outputvol'
    )


    smooth_only = traits.Bool(
        argsstr='--smooth-only',
        desc='smooth and save, do not compute fwhm (--so)'
    )

    mask_file = File(
        argstr='--mask %s',
        desc='binary mask'
    )

    mask_thresh = traits.Float(
        argstr='--mask-thresh %f',
        desc='absolute threshold for mask (default is .5)'
    )

    auto_mask_thresh = traits.Float(
        argstr='--auto-mask %f',
        desc='threshold for auto mask'
    )


    nerode = traits.Int(
        argstr='--nerode %i',
        desc='erode mask n times prior to computing fwhm'
    )

    mask_inv = traits.Bool(
        argstr='--mask-inv',
        desc='invert mask'
    )
    
    out_mask = File(
        argstr='--out-mask %s',
        desc='save final mask'
    )

    x = File(
        argstr='--X %s',
        desc='matlab4 detrending matrix'
    )
   
    detrend = traits.Int(
        argstr='--detrend %i',
        desc='polynomial detrending (default 0)'
    )

    sqr = traits.Bool(
       argstr='--sqr',
       desc='compute square of input before smoothing'
    )

    fwhm = traits.Float(
            argstr='--fwhm %f',
            desc='smooth BY fwhm before measuring'
        )
   
    gstd = traits.Float(
        argstr='--gstd %f',
        desc='same as --fwhm but specified as the stddev'
    )
   
    median = traits.Int(
        argstr='--median %i',
        desc='perform median filtering instead of gaussian'
    )

    to_fwhm = traits.Float(
        argstr='--to-fwhm %f',
        desc='smooth TO fwhm'
    )

    to_fwhm_tol = traits.Float(
        argstr='--to-fwhm-tol %f',
        desc='smooth to fwhm +/- tol (def .5mm)'
    
    )
    
    to_fwhm_nmax = traits.Int(
        argstr='--to-fwhm-nmax %i',
        desc='maximum number of iterations (def 20)'
    )
    
    to_fwhm_file = File(
        argstr='--to-fwhm-file %s',
        desc='save to-fwhm params in file'
    )
    
    sum_file = File(
        argstr='--sum %s',
        desc='summary/log'
    )
    dat = File(
        argstr='--dat %s',
        desc='only the final fwhm estimate'
    )
    
    synth = traits.Bool(
        argstr='--synth',
        desc='Create synthetic data (???)'
    )
    
    synth_frames = traits.Int(
        argstr='--synth-frames %i',
        desc='number of synthetic frames, default is 10'
    )
    
    
    nframesin = traits.Int(
        argstr='--nframesmin %i',
        desc='require at least this many frames'
    )
    
    ispm = traits.Bool(
        argstr='--ispm',
        desc='input is spm-analyze. Set --i to stem.'
    )
    in_nspmzeropad = traits.Int(
        argstr='--in_nspmzeropad %i',
        desc='size of zero-padding for spm-analyze'
    )
    
    nthreads = traits.Int(
        argstr='--nthreads %i',
        desc='Set OPEN MP threads'
    )

class SmoothVolOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Smoothed image")
    

class SmoothVol(FSCommand):
    
    """Smooth volume using mri_fwhm.

    Examples
    --------
    >>> smoothvol = SmoothVol()
    >>> smoothvol.inputs.in_file = 'img.nii'
    >>> smoothvol.inputs.out_file = 'smoothed_img.nii'
    >>> smoothvol.inputs.fwhm = 5
    >>> smoothvol.cmdline == 'mri_fwhm --fwhm 5.000000 --i img.nii --o smoothed_img.nii'
    """

    _cmd = "mri_fwhm"
    input_spec = SmoothVolInputSpec
    output_spec = SmoothVolOutputSpec

    def _format_arg(self, name, spec, value):       
        return super(SmoothVol, self)._format_arg(name, spec, value)
    

class GTMSegInputSpec(FSTraitedSpec):
    
    subject_id = traits.String(
        argstr="--s %s",
        desc="subject id",
        mandatory=True
    )

    xcerseg = traits.Bool(
        argstr="--xcerseg",
        desc="run xcerebralseg on this subject to create apas+head.mgz"
    )

    out_file = File(
        argstr="--o %s",
        desc="output volume relative to subject/mri (default is gtmseg.mgz)"
    )

    usf = traits.Int(
        argstr="--usf %i",
        desc="upsampling factor (default is 2)"
    )

    subsegwm = traits.Bool(
        argstr="--subsegwm",
        desc="subsegment WM into lobes (default)"
    )
    
    keep_hypo = traits.Bool(
        argstr="--keep-hypo",
        desc="do not relabel hypointensities as WM when subsegmenting WM"
    )
    
    keep_cc = traits.Bool(
        argstr="--keep-cc",
        desc="do not relabel corpus callosum as WM"
    )

    dmax = traits.Float(
            argstr="--dmax %f",
            desc="distance threshold to use when subsegmenting WM (default is 5)"
        )

    ctx_annot = traits.Tuple(
        traits.String,
        traits.Int,
        traits.Int,
        argstr="--ctx-annot %s %i %i",
        desc="annot lhbase rhbase : annotation to use for cortical segmentation (default is aparc 1000 2000)"
    )

    wm_annot = traits.Tuple(
        traits.String,
        traits.Int,
        traits.Int,
        argstr="--wm-annot %s %i %i",
        desc="annot lhbase rhbase : annotation to use for WM segmentation (with --subsegwm, default is lobes 3200 4200)"
    )

    output_usf = traits.Int(
        argstr="--output-usf %i",
        desc="set output USF different than USF, mostly for debugging"
    )

    head = traits.String(
        argstr="--head %s",
        desc="use headseg instead of apas+head.mgz"
    )

    subseg_cblum_wm = traits.Bool(
        argstr="--subseg-cblum-wm",
        desc="subsegment cerebellum WM into core and gyri"
    )     

    no_pons = traits.Bool(
        argstr="--no-pons",
        desc="do not add pons segmentation when doing ---xcerseg"
    )
    
    no_vermis = traits.Bool(
        argstr="--no-vermis",
        desc="do not add vermis segmentation when doing ---xcerseg"
    )

    ctab = File(
        exists=True,
        argstr="--ctab %s",
        desc="colortable"
    )
    no_seg_stats = traits.Bool(
        argstr="--no-seg-stats",
        desc="do not compute segmentation stats"
    )    


class GTMSegOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="GTM segmentation")


class GTMSeg(FSCommand):
    """create an anatomical segmentation for the geometric transfer matrix (GTM).

    Examples
    --------
    >>> gtmseg = GTMSeg()
    >>> gtmseg.inputs.out_file = 'gtmseg.nii'
    >>> gtmseg.inputs.subject_id = 'subjec_id'
    >>> gtmseg.cmdline == 'gtmseg --o gtmseg.nii --s subject_id'

    """

    _cmd = "gtmseg"
    input_spec = GTMSegInputSpec
    output_spec = GTMSegOutputSpec

    def _format_arg(self, name, spec, value):       
        return super(GTMSeg, self)._format_arg(name, spec, value)


class MRTMInputSpec(GLMFitInputSpec):
    
    mrtm = traits.Tuple(
        traits.String,
        traits.String,
        mandatory=True,
        argstr="--mrtm1 %s %s",
        desc="RefTac TimeSec : perform MRTM1 kinetic modeling"
    )


class MRTMOutputSpec(GLMFitInputSpec):
    r2p = File(desc="estimate of r2p parameter")


class MRTM(GLMFit):
    """Perform MRTM1 kinetic modeling.

    Examples
    --------
    >>> mrtm = MRTM()
    >>> mrtm.inputs.in_file = 'tacs.nii'
    >>> gtmseg.inputs.mrtm = ('ref_tac.dat', 'timing.dat')
    >>> mrtm.inputs.glmdir = 'mrtm'
    >>> mrtm.cmdline == 'mri_glmfit --glmdir mrtm --y tacs.nii --mrtm1 ref_tac.dat timing.dat'
    """

    _cmd = "mri_glmfit"
    input_spec = MRTMInputSpec
    output_spec = MRTMOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        return outputs


class MRTM2InputSpec(GLMFitInputSpec):
    
    mrtm2 = traits.Tuple(
        traits.String,
        traits.String,
        traits.Float,
        mandatory=True,
        argstr="--mrtm2 %s %s %f",
        desc="RefTac TimeSec k2prime : perform MRTM2 kinetic modeling"
    )
    
    _ext_xor = ['nii', 'nii_gz']
    nii = traits.Bool(
        argstr='--nii',
        desc='save outputs as nii',
        xor=_ext_xor
    )
    nii_gz = traits.Bool(
        argstr='--nii.gz',
        desc='save outputs as nii.gz',
        xor=_ext_xor
    )


class MRTM2OutputSpec(GLMFitInputSpec):
    bp = File(desc="BP estimates")


class MRTM2(GLMFit):
    """Perform MRTM2 kinetic modeling.

    Examples
    --------
    >>> mrtm = MRTM()
    >>> mrtm.inputs.in_file = 'tacs.nii'
    >>> gtmseg.inputs.mrtm = ('ref_tac.dat', 'timing.dat', 'k2prime.dat')
    >>> mrtm.inputs.glmdir = 'mrtm2'
    >>> mrtm2.cmdline == 'mri_glmfit --glmdir mrtm2 --y tacs.nii --mrtm2 ref_tac.dat timing.dat k2prime.dat'
    """

    _cmd = "mri_glmfit"
    input_spec = MRTM2InputSpec
    output_spec = MRTM2OutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.nii_gz):
            ext = '.nii.gz'
        if isdefined(self.inputs.nii):
            ext = '.nii'
        else:
            ext = '.mgh'            
        outputs['bp'] = join(self.inputs.glm_dir, 'bp' + ext)
        return outputs


class GCAMInputSpec(FSTraitedSpec):
    gcam = traits.Tuple(
            File,
            traits.String,  # can be either file or 0, not sure how to implement
            traits.String,
            traits.String,
            traits.String,
            traits.Int(min=0, max=5),
            traits.String,
            argstr="--gcam %s %s %s %s %s %i %s",
            desc="mov srclta gcam dstlta vsm interp out \
             srclta, gcam, or vsm can be set to 0 to indicate identity \
             direction is automatically determined from srclta and dstlta \
             interp 0=nearest, 1=trilin, 5=cubicbspline"
        )

class GCAMOutputSpec(TraitedSpec):
    out_file = File('input in target space')

class GCAM(FSCommand):
    """Apply a transform computed by mri_cvs_register.

    Examples
    --------
    >>> gcam = GCAM()
    >>> gcam.inputs.gcam = ('mov.nii', 'src.lta', 'gcam', 'dst.lta', 0, 5, out.nii'
    >>> gcam.cmdline == 'mri_vol2vol --gcam mov.nii src.lta gcam dst.lta 0 5 out.nii'
    """
    
    _cmd = "mri_vol2vol"
    input_spec = GCAMInputSpec
    output_spec = GCAMOutputSpec

#%% Specific processing

def get_timing(json_file):
    
    """
    Extract timing of dynamic PET data from JSON file 
    
    Arguments
    ---------
    json_file: string
        path to JSON file
    fout: string
        path to output average volume
    
    Return
    ---------
    frames_start: numpy array of float
        list of start time for each frame (sec) 
    frames_duration: numpy array of int
        list of frame duration (sec)
        
    """

    with open(json_file, 'r') as f:
        info = json.load(f)
    frames_start = np.array(info['FrameTimesStart'], dtype=float)
    frames_duration = np.array(info['FrameDuration'], dtype=float)
    return frames_start, frames_duration    

def create_weighted_average_pet(fin, json_file, fout, frames=None):
    
    """
    Create a time-weighted average of dynamic PET data using mid-frames
    
    Arguments
    ---------
    fin: string
        path to input dynamic PET volume
    fout: string
        path to output average volume
    frames: list of integers
        list of frames to be used for computing the average (indices are 0-based)
    """     
      
    if not isfile(fout):
        img = nib.load(fin)        
        data = img.get_fdata()

        frames_start, frames_duration = get_timing(json_file)
        
        # Check that the specified frame interval, if any, is valid
        if frames is None:
            frames = range(data.shape[-1])
        else:
            if frames[0] < 0:
                raise ValueError('The indice of of first frame needs to be equal or larger than 0')
            if frames[-1] >= data.shape[-1]:
                raise ValueError('The indice of of last frame needs to less than %i' % data.shape[-1])

        mid_frames = frames_start + frames_duration/2
        wavg = np.trapz(data[..., frames], dx=np.diff(mid_frames[frames]), axis=3)/np.sum(mid_frames)
        print('Saving average to ' + fout)
        nib.save(nib.Nifti1Image(wavg, img.affine), fout)
    else:
        print('File ' + fout + ' already exists. Skipping.')
        

def freeview_QC(cmd, png_concat, viewports=['sagittal', 'coronal', 'axial']):
    
    """
    Create images for quality control using freeview
    
    Arguments
    ---------
    cmd: string
        command to be executed by freeview
    png_concat: string
        path to file concatenating all views
    viewports: list of string
        views to be visualized (saggital, coronal, axial)
    """  
    
    # Handle strings
    if not isinstance(viewports, list):
        viewports = [viewports]
    
    pngs = []
    for viewport in viewports:
        png_out = join('/tmp', viewport + '.png')  # might conflict with parallel execution
        pngs += [png_out]
        cmd = ' '.join([cmd,
                '-viewport', viewport,
                '-ss ' + png_out,
                '-quit'
            ])
        run(cmd)
    
    # Concat images        
    cmd = ' '.join(['convert', ' '.join(pngs), '-append', png_concat])
    run(cmd)
    list(map(os.remove, pngs))  
    
def visual_coreg_QC(anat, pet, png, color_range=[50, 800], opacity=0.3):
    
    """
    Create images for quality control of the aligment between an average
    PET image and an anatomical image
    
    Arguments
    ---------
    anat: string
        path to anatomical image
    pet: string
        path to PET image
    png: string
        path to output image
    """  
    
    # Check if anat image exists
    if not isfile(anat):
        raise ValueError('Anatomical image does not exist. ' + anat)
        return

    # Check if reference image exists
    if not isfile(pet):
        raise ValueError('PET image does not exist. ' + pet)
            
    if not isfile(png):
        cmd = ' '.join([
                'freeview', anat,
                pet + ':colorscale=%i,%i:colormap=jet:opacity=%f' % (
                    color_range[0], color_range[1], opacity
                )
            ])
        freeview_QC(cmd, png)
    else:
        print('QC Image already exists. Skipping. ' + png)


def visual_gtmseg_QC(recon_dir, png):
    
    """
    Create images for quality control of the gtmseg segmentation
    
    Arguments
    ---------
    recon_dir: string
        path to the subject's recon directory
    png: string
        path to output image
    """  
    
    gtmseg = join(recon_dir, 'mri', 'gtmseg.mgz')
    
    if isfile(gtmseg) and not isfile(png):        
        anat = join(recon_dir, 'mri', 'norm.mgz')
        cmd = ' '.join([
                'freeview', anat,
                gtmseg + ':colormap=lut:opacity=0.3'
            ])
        freeview_QC(cmd, png)


def tac_file_exist(tacs_dir, labels_dct):
    
    """
    Check wether all output TAC files have been created
    
    Arguments
    ---------
    tac_dir: string
        path to directory containing TAC files
    labels_dct: dictionary
        dictionary specifying TAC files to be created
    """ 
    
    keys = list(labels_dct.keys())
    ext = [labels_dct[key]['ext'] for key in keys]
    return np.all([isfile(join(tacs_dir, k + '.' + e))
                   for k, e in zip(keys, ext)])


def save_np_array_to_fs(X, fout):
    
    """
    Save surface data from numpy array to to Nifti
    
    Arguments
    ---------
    X: numpy array
        data array to be save
    fout: string
        path to nifti file
    """  
    
    # Assert input shape
    if X.ndim > 2:
        raise ValueError('Input array has more than 2 dimensions')
    if X.ndim == 1:
        X = X.reshape([-1, 1])
    
    nib.save(nib.Nifti1Image(
                X.reshape([X.shape[0], 1, 1, X.shape[-1]]), np.eye(4)),
            fout)


def extract_vol_tacs(pet_file, out_dir, labels_dct):

    """
    Extract TACs from volume as specified in labels dictionary
    
    Arguments
    ---------
    pet_file: string
        path to input dynamic PET volume
    out_dir: string
        path to output directory
    labels_dct: dictionary
        dictionary specifying TAC files to be created
    """     
    
    assert_dir(out_dir)

    print('Loading data: ' + pet_file)
    data = nib.load(pet_file).get_fdata()
                
    for key in labels_dct.keys():
        
        print('Processing ' + key)
        
        ext = labels_dct[key]['ext']        
        fout_data = join(out_dir, key + '.' + ext)
        if not isfile(fout_data):

            labels = np.squeeze(nib.load(labels_dct[key]['file']).get_fdata())
            ids = labels_dct[key]['ids']

            # Extract mean tacs
            mean_tacs = []
            for label in ids:
                mask = np.isin(labels, label)
                mean_tacs += [np.mean(data[mask, :], axis=0)]
            mean_tacs = np.vstack(mean_tacs)
                
            
            print(fout_data)
            if ext == 'nii.gz':
                save_np_array_to_fs(mean_tacs, fout_data)
            elif ext == 'dat':
                np.savetxt(fout_data, mean_tacs.T, fmt='%0.8f')
            else:
                raise ValueError('Invalid extension ' + ext)
            

def extract_surf_tacs(surf_tacs, hemi, annot, fout):

    """
    Extract TACs from surface as specified FreeSurfer parcellations
    
    Arguments
    ---------
    surf_tacs: string
        path to input dynamic PET data sampled onto surface
    subject_id: string
        subject name in the recon directory (i.e., SUBJECTS_DIR)
    annot: string
        path to annot file for surface segmentation
    fout: string
        path to output file
    """     
        
    if not isfile(fout):
        data = nib.load(surf_tacs).get_fdata()
        labels, ctab, names = zip(nib.freesurfer.read_annot(annot))
        ids = np.unique(labels)
        ids = ids[ids != -1]
        
        # Extract mean tacs
        mean_tacs = []
        for label in ids:
            mask = np.isin(labels, label).reshape(-1)
            mean_tacs += [np.mean(data[mask, ...], axis=0).reshape(-1)]
        mean_tacs = np.vstack(mean_tacs)
    
        print(fout)        
        save_np_array_to_fs(mean_tacs, fout)


def create_mid_frame_dat(json_file, fout):

    """
    Extract timing of mid frames of dynamic PET data and save as text file for
    usage by mri_glmfit
    
    Arguments
    ---------
    json_file: string
        path to BIDS json PET file
    fout: string
        path to output file
    """  

    if not isfile(fout):
        frames_start, frames_duration = get_timing(json_file)
        mid_frames = frames_start + frames_duration/2    
        np.savetxt(fout, mid_frames, fmt='%0.1f')
    else:
        print('Midframe timing file already exists. Skipping.')
