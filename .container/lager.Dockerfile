#syntax=docker/dockerfile:1.2

## ========================================================================================
## STAGE1: EIC-dev base image
## ========================================================================================
FROM  eicweb/eic_dev:24.03-stable-default
LABEL maintainer="Sylvester Joosten <sjoosten@anl.gov>" \
      name="debian_base" \
      march="amd64"

## installing other dependencies (some of these should be folded into spack)
RUN  echo "  -------------------------------------------------"           \
  && echo "  1. Installing project dependency: hepmc2"                    \
  && echo "  -------------------------------------------------"           \
  && cd /tmp                                                              \
  && git clone --depth=1 https://eicweb.phy.anl.gov/monte_carlo/hepmc.git \
  && cmake -B build -S hepmc -Dmomentum:STRING=GEV -Dlength:STRING=CM  \
  && cmake --build build -j12 -- install                                  \
  && rm -rf build hepmc                                                   \
  && echo "  -------------------------------------------------"           \
  && echo "  2. Installing project dependency: PHOTOS++"                  \
  && echo "  -------------------------------------------------"           \
  && cd /tmp                                                              \
  && git clone --depth=1 https://eicweb.phy.anl.gov/monte_carlo/photospp.git \
  && cmake -B build -S photospp                                           \
  && cmake --build build -j12 -- install                                  \
  && rm -rf build photospp                                                \
  && echo "  -------------------------------------------------"           \
  && echo "  3. Installing: jpacPhoto"                                    \
  && echo "  -------------------------------------------------"           \
  && cd /tmp                                                              \
  && git clone --depth=1 -b v1.0.1                                        \
         https://eicweb.phy.anl.gov/monte_carlo/jpacPhoto.git jpac        \
  && cmake -B build -S jpac                                               \
  && cmake --build build -j12 -- install                                  \
  && rm -rf build jpac

## finally install lAger in separate layer
ARG REPO_URL="https://eicweb.phy.anl.gov/monte_carlo/lager.git"
ARG LAGER_VERSION=master
RUN  echo "  -------------------------------------------------"           \
  && echo "  FINALLY: Installing: LAGER"                                  \
  && echo "  -------------------------------------------------"           \
  && cd /tmp                                                              \
  && git clone -b $LAGER_VERSION --depth=1 $REPO_URL lager                \
  && cmake -B build -S lager                                              \
  && cmake --build build -j12 -- install                                  \
  && pushd lager                                                          \
  && echo " - LAGER: ${LAGER_VERSION}-$(git rev-parse HEAD)"              \
           >> /etc/lager_info                                             \
  && popd                                                                 \
  && rm -rf build lager

## extra postinstall hooks
RUN echo " - lager: nightly-jug_dev-based" >> /etc/jug_info               \
 && ls /etc/profile.d/ && rm -f /etc/profile.d/debuginfod.*
