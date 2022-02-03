# PETPipeline
Repository to showcase a pipeline for pre-processing PET data using Nipype workflows
<!-- Table of Contents --> 
<h2 id="table-of-contents"> :open_book: Table of Contents </h2>
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project"> About The Project</a></li>
    <li><a href="#prerequisites"> Prerequisites</a></li>
    <li><a href="#repository-overview"> Repository Overview</a></li>
    <li><a href="#installation">  Installation</a></li>
    <li><a href="#usage"> Usage</a></li>
    <li><a href="#configuration"> Configuration</a></li>
    <li><a href="#output-structure"> Output Structure</a></li>
    <li><a href="#license"> License</a></li>
  </ol>
</details>

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)

<!-- About the Project -->
<h2 id="about-the-project"> :memo: About the Project </h2>

*Positron Emision Tomography* (PET) is an advanced imaging technique used in the field of neuroimaging to quantify the concentration of molecular targets in the brain using intravenously injected radio-active tracers. The images acquired through PET undergo a series of preprocessing steps before they can be used for further analysis. Previous investigations have shown that differences in preprocessing strategies, study design, and statistical analysis among other variables may result in different conclusions.
This inturn challenges the reproducibility of the results and derived neurobiological interpretations. 
This repository showcases a robust and automated reproducible workflow for pre-processing of PET volumes built around Nipype, that can be applied across data from various sources to minimize unexplained variance originating from the differences in methodology. The pipeline includes various preprocessing steps commonly used for preprocessing of PET images such as 
*Motion Correction,* 
*Co-Registration*, 
*Delineation of Volumes of Interest*, 
*Partial Volume Correction* and 
*Pharmacokinetic Modelling.* 

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)


<!-- Prerequisites -->
<h2 id = "prerequisites"> :fork_and_knife: Prerequisites </h2>

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

The following open source packages are used in this project. For further information on them, check out their linked github repositories. 

* [Nipype](https://github.com/nipy/nipype) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/nipype)
* [PyBids](https://github.com/bids-standard/pybids) [![PyPI version](https://badge.fury.io/py/pybids.svg)](https://badge.fury.io/py/pybids)
* [Nibabel](https://nipy.org/nibabel/) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/nibabel)

The following neuroimaging software have been used for the development of this project:

* [Freesurfer](https://surfer.nmr.mgh.harvard.edu/)
* [PetSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/PetSurfer) 
* [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) *(FMRIB Software Library)*  

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)

<!-- Repository Overiview -->
<h2 id="repository-overview"> :world_map: Repository Overview </h2>

```
PETPIPELINE
├── LICENSE
├── petpipeline
│   ├── config.py
│   ├── config.yaml
│   ├── __init__.py
│   ├── main.py
│   ├── PETPipeline.py
│   ├── README
│   └── utils.py
├── README.md
├── setup.cfg
└── setup.py
```

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)

<!-- Installation -->
<h2 id="installation"> :gear: Installation </h2>

<h3 id="getting-started">Getting Started: Preparing the system</h3>

In order to run the pipeline, you can choose to either set up a docker container in the interactive mode or set up the requirements locally on your machine. 
Depending on the option you choose, you can follow along the section <a href="#environment-prep"> Preparing the Environment</a> section onwards to ensure that you have 
all the prerequisites needed to run the pipeline. 

### Setting up a Docker container

[Docker](https://www.docker.com/) is an open-source project for automating the deployment of applications as portable, self-sufficient containers that can run on the cloud. It provides a platform as a service products that use OS-level virtualization to deliver software in packages called containers. Containers are isolated from one another and bundle their own software, libraries and configuration files. They utilize the resources provided by the host machine. 
To get docker up and running on your system, you can follow along the instructions on the [webpage](https://docs.docker.com/engine/install/).

If you have docker set up on your system, then on a terminal run the following command to set up and run a container in an interactive mode. You can also attach the data you want to run the pipeline on, to your container using the `-v` flag as shown. This container uses an `ubuntu` image. :

```
docker container run  -v <path_to_data_in_host_machine>:<path_to_data_in_container>  --name <name_of_container> -it ubuntu
```

This will run a new container with the specified name and allow you to execute commands on a terminal inside the container. 
                                         
### Local Installation
Ensure that you have Linux or Mac OS up and running. Windows users can use [Virtual Box](https://www.virtualbox.org/wiki/Downloads) to create a [Virtual Machine](https://brb.nci.nih.gov/seqtools/installUbuntu.html) running a Linux OS.


<h2 id="environment-prep"> Preparing the Environment</h3>
To prepare the environment and install the required packages, run the commands below. 

Update and upgrade the environment using:

```
apt-get update && apt-get upgrade
```

Install `curl` and `wget` for downloading packages online. 

```
apt install curl
apt install wget
```

Install git 
```
apt install git 
```

Install a text editor: VIM, or use one if already installed.
```
apt install vim 
```

<h3 id="anaconda"> Install Anaconda</h3>

[Anaconda](https://www.anaconda.com/) is a distribution of the Python and R programming languages for scientific computing (data science, machine learning applications, large-scale data processing, predictive analytics, etc.), that aims to simplify package management and deployment. The distribution includes packages suitable for Windows, Linux, and macOS.

We'll be using Anaconda to install the latest version of Python version and its packages. 

To download anaconda, in the terminal, run the following command:

```
curl -O https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

To install, run:
```
bash Anaconda3-2021.11-Linux-x86_64.sh
```

Follow along the instructions as displayed on the terminal to successfully complete the installation. 

To verify if python as been installed successfully, on a terminal run `python`. 
It should run python and display the version of Python installed therby verifying the installation. 

<h3 id="neurodebian"> Install Neurodebian</h3>

[Neurodebian](http://neuro.debian.net/#get-neurodebian) provides a large collection of popular neuroscience research software for the Debian operating system as well as Ubuntu and other derivatives.
It makes the installation process of the neuroimaging software such as FSL easy and hassle-free. 

To download and install neurodebian, go to the [Neurodebian](http://neuro.debian.net/#get-neurodebian) webpage. 
Select the appropriate operating system and download server, and choose the `All Software` option. 
This will then generate the necessary commands needed for installation as shown below. On a terminal run these commands to complete the installation process.
For testing purposes, we've used Ubuntu 20.4 "Focal Fossa".

```
wget -O http://neuro.debian.net/lists/focal.de-fzj.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9
```

Update the package index using:
```
apt-get update
```

Now we can go ahead and install the neuroimaging software required!

<h3 id="neuroimaging-software"> Install Neuroimaging software </h3>

### **FSL**

[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/) *(FMRIB Software Library)*  is a comprehensive library of analysis tools for FMRI, MRI and DTI brain imaging data. 
Various tools that FSL provides are described in the [overview](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslOverview). 

We can install FSL using neurodebian using the following command:

```
apt-get install fsl
```

Edit the `bashrc` file using a text-editor:

```
vi ~/.bashrc
```

Add the following lines to it. This will add the path to `fsl` to the `PATH` environment variable:

```
#FSL
FSLDIR=/usr/share/fsl
. ${FSLDIR}/5.0/etc/fslconf/fsl.sh
PATH=${FSLDIR}/5.0/bin:${PATH}
export FSLDIR PATH
```

Source the `~/.bashrc` by running the following command on the terminal:

```
source ~/.bashrc
```

In order to verify that the installation has been completed successfully, on a terminal run the command `fsl`. 
It should display a list of options showing how to run `fsl`, thereby verifying the installation. 

### **Freesurfer**

[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) is an open source software suite for processing and analyzing (human) brain MRI images.
To download and install freesurfer, visit the [downloads](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) page and copy the link address for the appropriate version depending on your operating system. Also, register for a license for freesurfer on this [link](https://surfer.nmr.mgh.harvard.edu/registration.html). This should generate a license.txt file. 

On a terminal, run the following command with the copied link address following the -O- flag to start the download process:

```
curl -O https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz
```

Then run the following command to unzip and install freesurfer: 
```
tar -C /usr/local -xzvf freesurfer-Linux-<platform>-<release>-full.tar.gz
```

Freesurfer would now be installed in /usr/local/freesurfer. You can copy the contents of the `license` file to this directory using a text editor such as `vim`. 


Finally, edit the `bashrc` using a text editor:

```
vi ~/.bashrc 
```

Add the following lines to it: 
```
export FREESURFER_HOME=/usr/local/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
```

Source `~/.bashrc` by running the following command on the terminal:
```
source `~/.bashrc`
```

<h3 id="nipype-petsurfer"> Install Nipype & PETSurfer interface for Nipype </h3>

Install the [PETSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/PetSurfer) interface to Nipype from the [add_pet_surfer](https://github.com/mnoergaard/nipype/tree/add_pet_freesurfer) branch from the following [repository](https://github.com/mnoergaard/nipype/tree/add_pet_freesurfer):

```
mkdir petsurfer
cd petsurfer
git clone https://github.com/mnoergaard/nipype.git --branch add_pet_freesurfer
cd nipype
pip install .
```

This should install Nipype as well as the PETSurfer interface to Nipype. 
To verify the installation run `python` in the terminal then import `nipype` and `petsurfer` using the following commands:

 ```
 import nipype
 from nipype.interfaces.freesurfer import petsurfer
 ```
If the installation was successful, you should be able to import them without errors. 

Now that we have all the requirements in place, we can go ahead and install the PetPipline!
 
<h3 id="install-pipeline"> Install Petpipeline </h3>

Clone the repository using 

```
git clone https://github.com/openneuropet/petpipeline.git
```

Install python packages using requirements.txt

This should install python packages `pybids` and `pyyaml`

``` 
pip install requirements.txt
```

Specify the required directories and configuration parameters for the different preprocessing steps in the `config.yaml` file as explained in the <a href="#configuration"> configuration section. </a>

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)


<!-- Usage -->
<h2 id="usage"> :rocket: Usage </h2>
  
In order to run the pipeline ensure that your dataset is compatible with [BIDS](https://bids.neuroimaging.io/). You can validate the dataset using the [BIDS Validator](https://bids-standard.github.io/bids-validator/). 

An outline of a sample PETCimbi dataset from [OpenNeuro](https://openneuro.org/) is shown below. It is a representation of the standard your dataset should comply with.
```
PET_CIMBI
├── CHANGES
├── dataset_description.json
├── participants.json
├── participants.tsv
├── README
├── sub-01
│   ├── ses-baseline
│   │   ├── anat
│   │   │   ├── sub-01_ses-baseline_T1w.json
│   │   │   └── sub-01_ses-baseline_T1w.nii
│   │   └── pet
│   │       ├── sub-01_ses-baseline_pet.json
│   │       └── sub-01_ses-baseline_pet.nii.gz
│   └── ses-rescan
│       ├── anat
│       │   ├── sub-01_ses-rescan_T1w.json
│       │   └── sub-01_ses-rescan_T1w.nii
│       └── pet
│           ├── sub-01_ses-rescan_pet.json
│           └── sub-01_ses-rescan_pet.nii.gz
```

Navigate to the petpipeline sub-directory 
```
cd petpipeline
```

Run the pipeline using:
```
python main.py --config config.yaml
```

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)

<!-- Configuration -->
<h2 id="configuration"> :hammer_and_wrench: Configuration </h2>

The `config.yaml` file allows you to configure the paths to the various directories as well as configure certain parameters for each of the preprocessing steps
within the pipeline. 

First, you would need to set the paths to the following directories:
```
environment:
  experiment_dir: '' # directory for running the experiments
  output_dir: 'derivatives/' # directory where you would want to store your outputs
  working_dir: 'working_dir/' # working directory (used by nipype) for the experiments
  data_dir: '' # path to the data directory for your experiments
 ```
 
 You can then choose to run the pipeline with default parameters or specify additional parameters for each of the preprocessing steps in a similar format. 

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)


<!-- Output Structure -->
<h2 id ="output-structure"> :open_file_folder: Output Structure </h2>

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)

<!--References -->
<h2 id ="output-structure"> :books: References </h2>

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)


<!-- License -->
<h2 id="license"> :scroll: License </h2>

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/openneuropet/petpipeline/blob/main/LICENSE)

![-----------------------------------------------------](https://github.com/andreasbm/readme/blob/master/assets/lines/aqua.png)
