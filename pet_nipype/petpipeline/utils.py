import os

def create_mid_frame_dat(json_file):
    
    import os
    from os.path import join, isfile
    import numpy as np
    import json
    from nipype.utils.filemanip import split_filename

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

    pth, fname, ext = split_filename(json_file)
    print(fname)
    fout = "{}_time.dat".format(fname)
    with open(json_file, 'r') as f:
        info = json.load(f)
    frames_start = np.array(info['FrameTimesStart'], dtype=float)
    frames_duration = np.array(info['FrameDuration'], dtype=float)
    mid_frames = frames_start + frames_duration/2    
    np.savetxt(fout, mid_frames, fmt='%0.1f')
    
    return os.path.abspath(fout)


def compute_average(in_file, out_file=None):

    """
        A function to compute tehe average over all time frames
        for a given pet volume
        
        Parameters
        ----------
        in_file : str 
            input file path (str) for pet volume
        out_file : str 
            output file path (str) computed average

        Returns
        -------
        out_file : str 
            output file path (str) computed average
    """

    import os
    import nibabel as nib
    import numpy as np
    from nipype.utils.filemanip import split_filename

    pet_brain = nib.load(in_file)
    pet_brain_img = pet_brain.get_fdata()
    avg = np.mean(pet_brain_img, axis=3)
    pet_brain_frame = nib.Nifti1Image(avg, pet_brain.affine)
        
    new_pth = os.getcwd()
    pth, fname, ext = split_filename(in_file)
    pet_brain_filename = "{}_mean.nii.gz".format(fname)
    pet_brain_frame.to_filename(pet_brain_filename)

    return os.path.abspath(pet_brain_filename)

def compute_weighted_average(in_file, json_file, out_file=None): 

    """
        A function to compute a time weighted average over
        the time frames for a given pet volume

        Parameters
        ----------
        in_file : str 
            input file path (str) for pet volume
        out_file : str 
            output file path (str) computed average

        Returns
        -------
        out_file : str 
            output file path (str) computed average

    """   
        
    import numpy as np
    import os 
    import nibabel as nib
    import json
    from nipype.utils.filemanip import split_filename


    img = nib.load(in_file)        
    data = np.float32(img.get_fdata())

    with open(json_file, 'r') as f:
        desc = json.load(f)
        frames = np.float32(np.array(desc['FrameDuration'], dtype=float))

    data = np.sum(np.float32(data * frames),axis=3) / np.sum(frames)
      
    img_ = nib.Nifti1Image(data, img.affine)
            
    new_pth = os.getcwd()
    pth, fname, ext = split_filename(in_file)
    out_file = "{}_twa.nii.gz".format(fname)
    img_.to_filename(out_file)
    return os.path.abspath(out_file)    


def combine_file_paths(time_file, ref_file):

    """
        A function to return a tuple as a list
        using the input files provided

        Parameters
        ----------
        time_file : str 
            input file path (str) for pet volume
        ref_file : str 
            output file path (str) computed average

        Returns
        -------
        [(time_file, ref_file)] :
            [(time_file path (str), ref_file path)] 

    """

    import os
    return [(os.path.abspath(ref_file),os.path.abspath(time_file))]

def extract_value_from_file(in_file):
    with open(in_file) as f:
        val = float(f.readlines()[0].strip())
    return val

def combine_(ref_file, time_file, k2p_file):

    """
        A function to return a tuple as a list
        using the input files provided

        Parameters
        ----------
        time_file : str 
            input file path (str) for pet volume
        ref_file : str 
            output file path (str) computed average

        Returns
        -------
        [(time_file, ref_file, k2p)] :
            [(time_file path (str), ref_file path, k2p)] 

    """

    import os
    from utils import extract_value_from_file
    k2p = extract_value_from_file(k2p_file)

    return [(os.path.abspath(ref_file),os.path.abspath(time_file), k2p)]

def listify(*args):

    """
        A function that creates a tuple 
        of all the arguments passed, for 
        inputs to InputMultiPath Spec

        Parameters
        ----------
        *args: a list of arguments

        Returns
        -------
        [(*args)] :
                a single-element list of a tuple containing 
                the arguments passed

    """
    return [tuple(args)]

def assert_dir(dir_path):
    """
        A function that creates a 
        directory if it does not exist

        Parameters
        ----------
        dir_path: path to directory
    """
    full_path = os.path.abspath(dir_path)
    if not os.path.isdir(full_path):
        os.makedirs(full_path)
