import os
from dataclasses import dataclass
from nipype.interfaces.fsl.preprocess import MCFLIRTInputSpec
from nipype.interfaces.freesurfer.registration import MRICoregInputSpec
from nipype.interfaces.freesurfer.preprocess import ReconAllInputSpec

@dataclass
class _EnvConfig:

    """
        A configuration class for the pet pipeline

        Attributes
        ----------
        data_dir : str
            Path to data directory
        experiment_dir: str
            Path to experiments directory
        working_dir : str 
            Path to the working directory
        output_dir : str
            Path to output directory

    """
    
    data_dir: str
    experiment_dir: str
    working_dir: str
    output_dir: str

class _MotionCorrectionConfig(MCFLIRTInputSpec):

    """
        A configuration class for motion correction,
        extends from Nipype MCFlirtInputSpec
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
 
class _CoregistrationConfig(MRICoregInputSpec):
    """
        A configuration class for coregistration,
        extends from Nipype MRICoregInputSpec
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

class _ReconAllConfig(ReconAllInputSpec):
    """
        A configuration class for recon-all,
        extends from Nipype ReconAllInputSpec
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

@dataclass 
class _PartialVolumeCorrectionConfig:
    psf: int
    default_seg_merge: bool
    auto_mask: tuple
    km_ref: list
    km_hb: list
    no_rescale: bool
    save_input: bool