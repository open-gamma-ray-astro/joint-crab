FROM continuumio/miniconda3:4.5.4
RUN apt-get update && apt-get install -y build-essential

# define work dir
WORKDIR /usr/src/app/

# copy content
COPY . /usr/src/app/

# install environment
ENV DOCKER_INSIDE "True"
RUN conda install -q -y pyyaml
RUN python env.py

# install joint-crab package
RUN /bin/bash -c "pip install ."

# add external notebook functionality
RUN echo 'alias notebook="jupyter notebook --no-browser --ip=0.0.0.0 --allow-root"' >> ~/.bashrc
EXPOSE 8888

# define entry
CMD [ "/bin/bash" ]
