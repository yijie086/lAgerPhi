FROM  eicweb.phy.anl.gov:4567/containers/image_recipes/root_base:latest

LABEL maintainer="Sylvester Joosten<sjoosten@anl.gov>" \
      name="liege"        \
      group="monte_carlo" \
      march="native"      \
      basedist="ubuntu"   \
      base="root_base"

ARG APP_VERSION

ENV PYTHONPATH=/usr/local/lib:$PYTHONPATH \
    LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

RUN echo ${APP_VERSION} && exit 1

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
  && echo "  2. Installing project dependency: PHOTOS++" \
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
  && echo "  3. Installing: LIEGE" \
  && echo "  -------------------------------------------------" \
  && cd /tmp \
  && git clone https://eicweb.phy.anl.gov/monte_carlo/liege.git \
  && git checkout ${APP_VERSION}
  && mkdir liege/build && cd liege/build \
  && cmake .. \
  && make -j20 \
  && make install \
  && cd /tmp \
  && rm -r /tmp/liege

RUN source /usr/local/bin/thisroot.sh \
  && echo "  -------------------------------------------------" \
  && echo "  4. Slimming down image" \
  && echo "  -------------------------------------------------" \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && echo "  -------------------------------------------------" \
  && echo "  Done building and installing software" \
  && echo "  -------------------------------------------------"
