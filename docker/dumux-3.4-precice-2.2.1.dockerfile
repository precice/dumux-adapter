FROM ubuntu:20.04

ARG DUNEVERSION=2.7
ARG DUMUXVERSION=3.4
ARG PRECICEVERSION=2.2.1
ARG PRECICEUBUNTU=focal
# Install DUNE and DuMuX in /home/dumux
WORKDIR	/home/dumux

ENV DEBIAN_FRONTEND noninteractive

# Prohibit /usr/doc files from being installed
ADD 01_nodoc /etc/dpkg/dpkg.cfg.d/01_nodoc

# Install dependencies from
RUN apt-get update
RUN apt-get install -y apt-utils vim htop git
RUN apt-get install -y libopenblas-dev libmetis-dev openmpi-bin libopenmpi-dev git gcc g++ wget cmake sudo libparmetis-dev libeigen3-dev libxml2-dev libboost-all-dev pkg-config

# Obtain DuMuX and Dune repositories
COPY checkout_repositories.sh .
RUN bash checkout_repositories.sh ${DUNEVERSION} ${DUMUXVERSION}

# Install DUNE and DUMUX
COPY cmake-test.opts .
# Make dunecontrol available
ENV PATH /home/dumux/dune-common/bin:${PATH}
RUN dunecontrol --opts=cmake-test.opts all

# Make DUNE and DUMUX accessible
ENV DUNE_CONTROL_PATH=/home/dumux/dune-common/dune.module:/home/dumux/dune-geometry/dune.module:/home/dumux/dune-grid/dune.module:/home/dumux/dune-localfunctions/dune.module:/home/dumux/dune-istl/dune.module:/home/dumux/dune-subgrid/dune.module:/home/dumux/dumux/dune.module

ENV PYTHONPATH /home/dumux/dumux/bin/testing:${PYTHONPATH}

#Install preCICE
RUN wget https://github.com/precice/precice/releases/download/v${PRECICEVERSION}/libprecice2_${PRECICEVERSION}_${PRECICEUBUNTU}.deb -O libprecice${PRECICEVERSION}-${PRECICEUBUNTU}.deb \
	&& apt install -y ./libprecice${PRECICEVERSION}-${PRECICEUBUNTU}.deb \
	&& rm -f ./libprecice${PRECICEVERSION}-${PRECICEUBUNTU}.deb

# Cleanup
RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/* && rm -rf /tmp/*

# Create user
ARG USER_NAME=dumux
ARG USER_HOME=/home/dumux
ARG USER_ID=1000
ARG GROUP_ID=1001

RUN groupadd -g $GROUP_ID "$USER_NAME"
RUN adduser \
    --home "$USER_HOME" \
    --uid $USER_ID \
    --gid $GROUP_ID \
    --disabled-password \
    "$USER_NAME"

RUN echo "$USER_NAME" ALL=\(root\) NOPASSWD:ALL > "/etc/sudoers.d/$USER_NAME" && \
    chmod 0440 "/etc/sudoers.d/$USER_NAME"

USER "$USER_NAME"
WORKDIR "$USER_HOME"

CMD ["bash"]

# Reset frontend for interactive use
ENV DEBIAN_FRONTEND=
