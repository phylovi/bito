FROM continuumio/anaconda3

RUN conda install -y -c anaconda \
        flex \
        gxx_linux-64 \
        make \
        pybind11 \
        scons

COPY . /libsbn
WORKDIR /libsbn
