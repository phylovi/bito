FROM matsengrp/conda-beagle

COPY environment.yml .

RUN /opt/conda/bin/conda env update -n libsbn -f environment.yml

WORKDIR /
ARG BEAGLE_PREFIX=/usr/local
