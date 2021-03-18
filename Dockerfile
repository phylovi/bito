ARG BASE=quay.io/matsengrp/conda-beagle:latest
FROM $BASE

COPY environment.yml .

RUN /opt/conda/bin/conda env create -f environment.yml

WORKDIR /
ENV BEAGLE_PREFIX /usr/local
ENV LD_LIBRARY_PATH /usr/local/lib
