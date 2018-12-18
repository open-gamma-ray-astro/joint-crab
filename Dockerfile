FROM continuumio/miniconda3:4.5.11
RUN apt-get update && apt-get install -y build-essential

# temporal folder to config environment
COPY binder/environment.yml tmp/
WORKDIR tmp/

RUN conda env create -f environment.yml
RUN echo 'conda activate joint-crab' >> ~/.bashrc

# define work dir
WORKDIR /usr/src/app/
ENV DOCKER_INSIDE "True"

# copy content
COPY . /usr/src/app/

# add external notebook functionality
RUN echo 'alias notebook="jupyter notebook --no-browser --ip=0.0.0.0 --allow-root"' >> ~/.bashrc
EXPOSE 8888

# define entry
CMD [ "/bin/bash" ]
