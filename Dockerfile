FROM continuumio/anaconda3

COPY . /libsbn
WORKDIR /libsbn

RUN conda env create -f environment.yml
RUN conda activate libsbn
RUN conda install -y gxx_linux-64
