FROM continuumio/anaconda3:2019.07

COPY environment.yml .

RUN /opt/conda/bin/conda create -n libsbn python=3.7 gxx_linux-64
RUN /opt/conda/bin/conda env update -n libsbn -f environment.yml
