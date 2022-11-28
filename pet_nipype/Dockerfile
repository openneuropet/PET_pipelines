FROM neurodebian:focal
ENV DEBIAN_FRONTEND=noninteractive
# update apt repos
RUN apt-get update && apt-get upgrade -y
# install global deps
RUN apt install curl wget git  fsl -y
RUN cat << EOF > ~/.bashrc \
    FSLDIR=/usr/share/fsl \
    . ${FSLDIR}/5.0/etc/fslconf/fsl.sh \
    PATH=${FSLDIR}/5.0/bin:${PATH} \
    export FSLDIR PATH \
    EOF \
    source ~/.bashrc

# collect and install anaconda
RUN curl -O https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh && bash Anaconda3-2021.11-Linux-x86_64.sh -b

# install freesurfer
RUN curl -O https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz
RUN cat << EOF > ~/.bashrc \
    export FREESURFER_HOME=/usr/local/freesurfer \
    source $FREESURFER_HOME/SetUpFreeSurfer.sh
RUN source ~/.bashrc

COPY . pet_nipype


