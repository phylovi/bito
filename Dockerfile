FROM continuumio/anaconda3

RUN conda install -y -c anaconda pybind11 scons make flex

COPY . /libsbn
WORKDIR /libsbn

#bison gxx_linux-64 numpy

# COPY . /bpp
# WORKDIR /bpp
# RUN python get_latest_bpp.py
# ENV LD_LIBRARY_PATH /usr/local/lib
