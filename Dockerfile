FROM continuumio/miniconda3

# set up conda and deps, per .travis.yml
RUN apt-get update -qq
RUN mkdir -p /usr/share/man/man1/  # required to prevent default-jre installation from crashing
RUN apt-get install -qq default-jre
RUN apt-get install -qq build-essential cmake

RUN conda config --set always_yes yes --set changeps1 no
RUN conda create -q -n py2 python=2.7

# copy source code from git
WORKDIR /mofid
COPY . /mofid

# compile openbabel, C++ analysis code, and Python package
RUN make init
RUN python set_paths.py
RUN pip install .
# Disabling conda environment for now and keeping with base Python
#RUN source activate py2
#RUN pip install .
#RUN conda deactivate

# then run the code, either interactively or via docker command line

