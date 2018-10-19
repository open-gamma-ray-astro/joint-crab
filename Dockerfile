FROM continuumio/miniconda3:4.5.11

RUN apt-get update && apt-get install -y build-essential

WORKDIR /usr/src/app/
COPY . /usr/src/app/

RUN conda env create -f config/environment.yml
RUN /bin/bash -c "source /opt/conda/bin/activate joint-crab"
RUN /bin/bash -c "python setup.py install"
RUN echo 'alias notebook="jupyter notebook --no-browser --ip=* --allow-root"' >> ~/.bashrc

ENV DOCKER_INSIDE "True"
EXPOSE 8888

CMD [ "/bin/bash" ]
