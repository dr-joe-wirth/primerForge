# syntax=docker/dockerfile:1

# Author: Joseph S. Wirth

# start with ubuntu base
FROM ubuntu:22.04

# install wget and python
RUN apt-get update
RUN apt-get -y install wget
RUN apt-get -y install python3.11 python3-pip python3-setuptools python3-dev

# install biopython
RUN pip3 install --upgrade pip
RUN pip3 install biopython

# install primerForge.py
RUN wget https://github.com/dr-joe-wirth/primerForge/archive/refs/tags/v0.2.0.tar.gz
RUN tar xzf v0.2.0.tar.gz
RUN rm v0.2.0.tar.gz
RUN mv /primerForge-0.2.0 /primerForge
RUN chmod +x /primerForge/primerForge.py

# add executables to the path
ENV PATH="${PATH}:/primerForge"

# create the volume and set it as the working directory
VOLUME ["/data"]
WORKDIR /data

# add an entry point so that primerForge.py can be called without entering the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
