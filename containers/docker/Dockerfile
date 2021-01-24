FROM  eicweb.phy.anl.gov:4567/containers/image_recipes/root_base:latest

LABEL maintainer="Sylvester Joosten<sjoosten@anl.gov>" \
      name="lager"        \
      group="monte_carlo" \
      march="native"      \
      basedist="ubuntu"   \
      base="root_base"

ARG APP_VERSION
ARG REPO_URL

ENV PYTHONPATH=/usr/local/lib:$PYTHONPATH \
    LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  1. Installing project dependency: hepmc" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone https://eicweb.phy.anl.gov/monte_carlo/hepmc.git \
  && mkdir build-hepmc && cd build-hepmc \
  && cmake ../hepmc -Dmomentum:STRING=GEV -Dlength:STRING=CM \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/*hepmc

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  2. Installing project dependency: hepmc3" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone https://gitlab.cern.ch/hepmc/HepMC3.git hepmc3 \
  && mkdir build-hepmc3 && cd build-hepmc3 \
  && cmake ../hepmc3 -DROOT_DIR=`root-config --prefix` \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/*hepmc3

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  3. Installing project dependency: PHOTOS++" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone https://eicweb.phy.anl.gov/monte_carlo/photospp.git \
  && mkdir photospp/build && cd photospp/build \
  && cmake .. \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/photospp

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  5. Installing: jpacPhoto" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone https://eicweb.phy.anl.gov/monte_carlo/jpacPhoto.git jpac \
  && cd jpac && git checkout v1.0.1 \
  && mkdir build && cd build \
  && cmake .. \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/jpac

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  5. Installing: LAGER" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone $REPO_URL lager \
  && cd lager \
  && git checkout $APP_VERSION \
  && mkdir build && cd build \
  && cmake .. \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/lager

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  6. Adding missing packages and slimming down image" \
  && echo "  -------------------------------------------------" \
  && apt-get -yqq update \
  && apt-get install -yqq parallel \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && echo "  -------------------------------------------------" \
  && echo "  Done building and installing software" \
  && echo "  -------------------------------------------------"
