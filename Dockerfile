# MOFid 2024
FROM python:3.11-bookworm

LABEL AUTHOR="darinweber@u.northwestern.edu"
LABEL version="2024"

WORKDIR /mofid
COPY . /mofid

RUN apt-get update -qq \ 
    && apt-get install -qq default-jre gcc-11 g++-11 cmake \
    && make init \
    && python set_paths.py \ 
    && pip install . \
