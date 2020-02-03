FROM  eicweb.phy.anl.gov:4567/containers/image_recipes/root_base:latest

LABEL maintainer="Sylvester Joosten<sjoosten@anl.gov>" \
      name="liege"        \
      group="monte_carlo" \
      march="native"      \
      basedist="ubuntu"   \
      base="root_base"

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
  && git clone https://eicweb.phy.anl.gov/monte_carlo/liege.git \
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


#

RUN source /usr/local/bin/thisroot.sh \
  && export PYTHONPATH=/usr/local/lib:$PYTHONPATH  \
  && git clone https://gitlab+deploy-token-1:FaNA-Yg4s7hpjvWPZnq8@eicweb.phy.anl.gov/upsilon/Pcsim.git \
  && mkdir Pcsim/build && cd Pcsim/build  \
  && cmake ../.  \
  && make -j20 \
  && make install


#-DCMAKE_CXX_FLAGS=" -march=haswell -O3 -mfma -malign-data=cacheline -finline-functions "
#&& wget -O- https://root.cern.ch/download/root_v6.14.06.source.tar.gz | tar -zxvf - \
#&& mv root-6.14.06 root_master \
#RUN which c++ && ls -lrth /usr/bin/c++ && cd /tmp/builds/root_build && make -j38 VERBOSE=1 && make install \
#      && cd /tmp && rm -rf /tmp/root_master && rm -rf /tmp/builds/root_build 

#RUN useradd -ms /bin/bash -d /opt/user user
#USER user
#WORKDIR /opt/bubble_user

##CMD ["-c" ]
#ENTRYPOINT ["/bin/bash"]

