# syntax=docker/dockerfile:1.4.0
FROM neurodebian:focal
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

# update apt repos
RUN apt-get update && apt-get upgrade -y
# install global deps
RUN apt install \
curl \
wget \
git \
fsl \
python3 \
python3-pip \
-y

# next we collect all of the other dependencies we need

# install freesurfer
RUN curl -O https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz && \
tar -C /usr/local -xzvf /freesurfer*

# now we configure our paths so our installs are reachable

# set up fsl path
RUN echo -e \ 
'FSLDIR=/usr/share/fsl\n. ${FSLDIR}/5.0/etc/fslconf/fsl.sh\nPATH=${FSLDIR}/5.0/bin:${PATH}\nexport FSLDIR PATH' \
>> ${HOME}/.bashrc && source ${HOME}/.bashrc

# Set up freesurfer on path
RUN echo -e \
"FREESURFER_HOME=/usr/local/freesurfer\nexport FREESURFER_HOME\nsource \$FREESURFER_HOME/SetUpFreeSurfer.sh\n" \
>> ${HOME}/.bashrc && source ${HOME}/.bashrc

# install nipype and petsurfer interface
RUN mkdir petsurfer && \
cd petsurfer && \
git clone https://github.com/mnoergaard/nipype.git --branch add_pet_freesurfer && \
cd nipype && \
pip install .

# Make processing directories
RUN mkdir /experiment_dir &&  \
mkdir /output_dir && \
mkdir /working_dir && \
mkdir data_dir

COPY . pet_nipype 
COPY petpipeline/config.yaml /config.yaml

RUN cd pet_nipype/ && pip install .

ENTRYPOINT [ "python3", "/pet_nipype/petpipeline/main.py" ]

CMD [ "--config", "/config.yaml", "--output_dir", "/output_dir", "--working_dir", "/working_dir", "--experiment_dir", "/experiment_dir", "--data_dir", "/data_dir" ] 