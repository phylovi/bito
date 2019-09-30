FROM continuumio/anaconda3

COPY environment.yml .
COPY beagle-recipe/ .

RUN /opt/conda/bin/conda env create -f environment.yml
RUN /opt/conda/bin/conda install -n libsbn -y gxx_linux-64
RUN /opt/conda/bin/conda build beagle-recipe
RUN /opt/conda/bin/conda install -n libsbn --use-local beagle-lib
