FROM continuumio/miniconda3:4.5.4
RUN apt-get update && apt-get install -y build-essential

# temporal folder
COPY env.py tmp/
COPY setup.py tmp/
COPY config/ tmp/config
COPY joint_crab/ tmp/joint_crab
WORKDIR tmp/

# install environment
RUN conda install -q -y pyyaml
RUN python env.py

# install joint-crab package
RUN /bin/bash -c "pip install ."

# add jointcrab user running the jupyter notebook process
ENV NB_USER jointcrab
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# copy content
COPY . ${HOME}

# setting ownerships
RUN chown -R ${NB_UID} ${HOME}

# start JupyterLab server in home dir
USER ${NB_USER}
WORKDIR ${HOME}
