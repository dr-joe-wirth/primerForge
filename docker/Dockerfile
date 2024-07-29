# syntax=docker/dockerfile:1

# Author: Joseph S. Wirth

# start with ubuntu base
FROM ubuntu:22.04

# install wget and python
RUN apt-get update && \
    apt-get -y install python3.11 python3-pip python3-setuptools python3.11-dev

# install conda
RUN wget -O conda_installer.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
RUN bash conda_installer.sh
RUN conda update conda

# install primerForge.py
RUN python3.11 -m pip install primerforge

# install isPcr
RUN conda install isPcr

# verify installation
RUN primerForge --help
RUN primerForge --check_install

# create the volume and set it as the working directory
VOLUME ["/data"]
WORKDIR /data

# add an entry point so that primerForge.py can be called without entering the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
