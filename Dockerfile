FROM continuumio/anaconda3

COPY . /libsbn
WORKDIR /libsbn

RUN conda install -y -c anaconda \
        flex \
        gxx_linux-64 \
        make \
        pybind11 \
        scons

RUN conda env create -f environment.yml
RUN conda activate libsbn
RUN conda install -y gxx_linux-64
