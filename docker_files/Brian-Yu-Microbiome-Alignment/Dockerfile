# Set the base image
FROM ubuntu:16.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER Brian Yu <brian.yu@czbiohub.org>

# Add packages, update image, and clear cache
RUN apt-get update
RUN apt-get install -y apt-utils 
RUN apt-get install -y build-essential wget zip unzip bzip2 git zlib1g-dev pkg-config make libbz2-dev python-pip

# RUN pip install --upgrade pip
# RUN pip install awscli --upgrade

# Install Anaconda3
# WORKDIR $HOME
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh -O /tmp/anaconda.sh
RUN /bin/bash /tmp/anaconda.sh -b -p ~/anaconda
RUN rm /tmp/anaconda.sh
RUN export PATH=~/anaconda/bin/conda:~/anaconda/bin/:$PATH
RUN ls
RUN /bin/bash ~/anaconda/bin/conda update -q -y conda
# Will install application specific conda environments in the aegea batch script

# End setup
# RUN apt-get clean
# CMD /bin/bash
# WORKDIR $HOME
