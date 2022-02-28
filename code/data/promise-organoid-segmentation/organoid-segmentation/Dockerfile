FROM nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu18.04

###################################################################################################
# SET UP SYSTEM
#
# Install OS-level packages
###################################################################################################
RUN apt update && apt upgrade -y && apt install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa -y && \
    apt update && \
    apt install -y \
    python3.6 \
    python3.6-dev \
    python3.6-venv \
    python3-distutils \
    ffmpeg \
    libsm6 \
    libxext6 \
    libxrender-dev \
    wget \
    gcc \
    gnupg

# Delete temp files to reduce image space
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

###################################################################################################
# CREATE VIRTUAL ENVIRONMENT
#
# For easier maintenance of files, all python libraries are installed into a virtual environment.
# This section also installs all libraries required by the AI to ensure the correct version is
# installed.
###################################################################################################
ENV VIRTUALENV=/.venv/myenv
RUN python3.6 -m venv $VIRTUALENV
ENV PATH=$VIRTUALENV/bin:$PATH
RUN pip install -U pip wheel
RUN pip install \
    tensorflow-gpu==2.6.2 \
    numpy==1.19.5 \
    scipy==1.5.4

###################################################################################################
# COPY FILES
###################################################################################################
WORKDIR /app
COPY weights/DNN_weights.h5 DNN_weights.h5
COPY image_preprocessor.py .
COPY main.py .
COPY entrypoint.sh .

ENTRYPOINT ["/bin/bash", "-c", "./entrypoint.sh"]