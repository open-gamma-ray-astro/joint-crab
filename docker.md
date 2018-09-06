## Docker

In order to use the `joint-crab` bundle in a virtual docker container you need to have [Docker](https://www.docker.com/community-edition) installed in your local machine.

### The `joint-crab` docker image

You can use a docker image that is stored in the [Gammapy DockerHub repository](https://hub.docker.com/u/gammapy/dashboard/). This docker image is atomatically downloaded when executing the `run` commands below. Alternatively you can build the docker image by your own with the [Dockerfile](Dockerfile) that we provide in this repository. In that case just type `docker build -t joint-crab .` at the top level of the joint-crab folder and remove `gammapy/` from the commands below.

You may reproduce the results runinng a docker container of the image in non-interactive mode or using an interactive bash session attached to the docker virtual environment.

### Non interactive mode

The following command will start up the docker container and reproduce the results of the paper using the docker virtual environment. The results are copied to your local filesystem in the folder `joint-crab-results` at your home directory.

    $ docker run --rm -it -v ~:/usr/src/app/host gammapy/joint-crab python /usr/src/app/make.py all


### Interactive bash session

    $ docker run --name joint-crab -it -p 9000:8888 -v ~:/usr/src/app/host gammapy/joint-crab

**Once the docker container running you land in the docker shell as a root user.**

We will use the prompt % to identify the commands in the docker container shell.

From this moment you are in the same executing environment that has been used to generate the results of this paper. You can now reproduce and inspect the results, execute the different notebooks or create new ones by re-use.

#### Reproduce the results

        % ./make.py all

The results are produced in folder `results` and copied to your local filesystem in the folder `joint-crab-results` at your home directory.

#### Execute notebooks

* Enter in the docker container shell

        % notebook

* Open a web browser and type

        http://localhost:9000

* Copy the token from the docker container shell into the form displayed by the browser and submit

#### Browse the docker container filesystem

* In a web browser

        http://localhost:9000

* In the docker container shell

        % cd results/figures
        % ls

#### Execute scripts or commands

* In the docker container shell

        % ./make.py --help   # See what is available

* [More details on how to reproduce the figures](analysis.md)

* Open an IPython session

         % ipython

#### Copy files from/into the docker container filesystem

* From the docker container to your home directory

        % cp -R notebooks host/

* From your home directory to the docker container

        % cp host/<myfile> .

*  In the browser

        http://localhost:9000

#### Finish your session and stop the docker container

        % exit

#### Recover a previous session re-starting the *join-crab* docker container

        $ docker start -ai joint-crab
