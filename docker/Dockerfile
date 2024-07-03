# syntax=docker/dockerfile:1

# Author: Joseph S. Wirth

# start with ubuntu base
FROM ubuntu:22.04

# install wget and python
RUN apt-get update && \
    apt-get -y install python3.11 python3-pip python3-setuptools python3.11-dev

# install primerForge.py
RUN python3.11 -m pip install primerforge

# create the volume and set it as the working directory
VOLUME ["/data"]
WORKDIR /data

# add an entry point so that primerForge.py can be called without entering the container
ENTRYPOINT ["/bin/bash", "-c", "exec \"$@\"", "-"]
