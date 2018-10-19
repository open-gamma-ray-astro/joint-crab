FROM continuumio/miniconda3:4.5.4
RUN apt-get update && apt-get install -y build-essential

# define work dir
WORKDIR /usr/src/app/

# copy content
COPY . /usr/src/app/

# install environment
RUN conda install -q -y pyyaml==3.12
RUN python env.py

# install joint-crab package
RUN /bin/bash -c "pip install ."
