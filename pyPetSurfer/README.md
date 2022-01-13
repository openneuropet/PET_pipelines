# Processing of PET data with FreeSurfer (PET Surfer)

This repository contains python code for the processing of PET data using the PET Surfer pipeline (https://surfer.nmr.mgh.harvard.edu/fswiki/PetSurfer).

The steps for performing regional, surface-based, and volume-based kinetic modeling using MRTM2 are outlined in **example.py**. This example presents the specific case of [11C]SB207145, and a test dataset for this tracer can be downloaded from [OpenNeuro](https://openneuro.org/datasets/ds001421). The FreeSurfer reconstruction and mask for the reference region are provided as derivatives in this dataset. The script will attempt to download the example dataset using [openneuro-cli](https://docs.openneuro.org/packages-openneuro-cli-readme), however, it can also be manually downloaded and extracted to a directory named ```ds001421-download``` located in the current directory from which the script is ran.

Required packages to run **example.py** are:
```
numpy
nibabel
nipype
matplotlib
```

## References

PET Surfer:
Greve DN, Svarer C, Fisher PM, et al. Cortical surface-based analysis reduces bias and variance in kinetic modeling of brain PET data. Neuroimage. 2014;92C:225-236. doi:10.1016/j.neuroimage.2013.12.021

MRTM:
Ichise M, Ballinger JR, Golan H, et al. Noninvasive quantification of dopamine D2 receptors with iodine-123-IBF SPECT. J Nucl Med. 1996;37(3):513-520.

MRTM2:
Ichise M, Liow J-S, Lu J-Q, et al. Linearized Reference Tissue Parametric Imaging Methods: Application to [11C]DASB Positron Emission Tomography Studies of the Serotonin Transporter in Human Brain. J Cereb Blood Flow Metab. 2003;23(9):1096-1112. doi:10.1097/01.WCB.0000085441.37552.CA
