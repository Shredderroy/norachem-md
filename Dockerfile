# ref https://github.com/tebeka/pythonwise/blob/master/docker-miniconda/Dockerfile
# FROM nvidia/cuda:11.0-devel-ubuntu18.04
FROM nvidia/cuda:11.0.3-devel-ubuntu20.04
# FROM nvidia/cuda:11.7.0-devel-ubuntu20.04


RUN bash
# System packages 
ENV DEBIAN_FRONTEND=nonintercative
RUN rm /etc/apt/sources.list.d/cuda.list 

RUN apt-get update --fix-missing 
RUN apt-get install -y curl
#RUN apt-get install -y vim

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV AMBERHOME=/miniconda
ENV PATH=/miniconda/bin:${PATH}

RUN conda update -y conda
COPY environment_docker.yml .
RUN conda env update --file environment_docker.yml --prune
#RUN conda env create -f environment_docker.yml 
RUN . /root/.bashrc
RUN conda init

COPY src /
COPY src/SIM.py .
COPY src/eval.sh .


ENTRYPOINT ["/eval.sh"]

