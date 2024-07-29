# syntax=docker/dockerfile:1

# Author: Joseph S. Wirth

# start with ubuntu base
FROM ubuntu:22.04

# install wget and python
RUN apt-get update && \
    apt-get -y install wget && \
    apt-get -y install python3.11 python3-pip python3-setuptools python3.11-dev

# install and configure conda
RUN wget -O /tmp/miniconda.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
RUN bash /tmp/miniconda.sh -b -p /opt/conda
ENV PATH /opt/conda/bin:$PATH
RUN rm /tmp/miniconda.sh
RUN conda update conda
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# install primerForge
RUN python3.11 -m pip install primerforge

# install isPcr
RUN conda install ispcr

# verify installation
RUN primerForge --help
RUN primerForge --check_install

# create the volume and set it as the working directory
VOLUME ["/data"]
WORKDIR /data

# add an entry point so that primerForge can be called without entering the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
