FROM continuumio/anaconda3

COPY environment.yml .

RUN /opt/conda/bin/conda env create -f environment.yml
RUN /opt/conda/bin/conda install -n libsbn -y gxx_linux-64
