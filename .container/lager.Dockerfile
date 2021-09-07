#syntax=docker/dockerfile:1.2

## ========================================================================================
## STAGE1: Debian base image
## ========================================================================================
FROM  amd64/debian:testing-20210408-slim as base
LABEL maintainer="Sylvester Joosten <sjoosten@anl.gov>" \
      name="debian_base" \
      march="amd64"

COPY bashrc /root/.bashrc

ENV CLICOLOR_FORCE=1                                                    \
    LANGUAGE=en_US.UTF-8                                                \
    LANG=en_US.UTF-8                                                    \
    LC_ALL=en_US.UTF-8

## Install additional packages. Remove the auto-cleanup functionality
## for docker, as we're using the new buildkit cache instead.
## We also install gitlab-runner, from the buster package (as bullseye is not available atm)
RUN --mount=type=cache,target=/var/cache/apt                            \
    rm -f /etc/apt/apt.conf.d/docker-clean                              \
 && apt-get -yqq update                                                 \
 && apt-get -yqq install --no-install-recommends                        \
        ca-certificates                                                 \
        curl                                                            \
        file                                                            \
        build-essential                                                 \
        g++-10                                                          \
        gcc-10                                                          \
        gdb                                                             \
        gfortran-10                                                     \
        git                                                             \
        gnupg2                                                          \
        iproute2                                                        \
        iputils-ping                                                    \
        locales                                                         \
        lua-posix                                                       \
        make                                                            \
        unzip                                                           \
        nano                                                            \
        vim-nox                                                         \
        less                                                            \
        clang-format                                                    \
        openssh-client                                                  \
        wget                                                            \
        ghostscript                                                     \
        gv                                                              \
        poppler-utils                                                   \
        parallel                                                        \
        time                                                            \
        valgrind                                                        \
 && localedef -i en_US -f UTF-8 en_US.UTF-8                             \
 && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100  \
 && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100  \
 && update-alternatives --install /usr/bin/gfortran gfortran            \
                                  /usr/bin/gfortran-10 100              \
 && cc --version                                                        \
 && curl -L                                                             \
    "https://packages.gitlab.com/install/repositories/runner/gitlab-runner/script.deb.sh" \
    | bash                                                              \
 && sed -i "s/bullseye/buster/"                                         \
           /etc/apt/sources.list.d/runner_gitlab-runner.list            \
 && apt-get -yqq update                                                 \
 && apt-get -yqq install --no-install-recommends                        \
        gitlab-runner                                                   \
 && apt-get -yqq autoremove                                             \
 && rm -rf /var/lib/apt/lists/*                                         

## ========================================================================================
## STAGE2: Builder image
## ========================================================================================

from base as builder

## instal some extra spack dependencies
RUN --mount=type=cache,target=/var/cache/apt                            \
    rm -f /etc/apt/apt.conf.d/docker-clean                              \
 && apt-get -yqq update                                                 \
 && apt-get -yqq install --no-install-recommends                        \
        python3                                                         \
        python3-pip                                                     \
        python3-setuptools                                              \
        tcl                                                             \
        uuid-dev                                                        \
        libfcgi-dev                                                     \
        x11proto-xext-dev                                               \
 && pip3 install boto3                                                  \ 
 && rm -rf /var/lib/apt/lists/*

## Setup spack
## parts:
ENV SPACK_ROOT=/opt/spack
ARG SPACK_VERSION="develop"
RUN echo "Part 1: regular spack install (as in containerize)"           \
 && git clone https://github.com/spack/spack.git /tmp/spack-staging     \
 && cd /tmp/spack-staging && git checkout $SPACK_VERSION && cd -        \
 && mkdir -p $SPACK_ROOT/opt/spack                                      \
 && cp -r /tmp/spack-staging/bin $SPACK_ROOT/bin                        \
 && cp -r /tmp/spack-staging/etc $SPACK_ROOT/etc                        \
 && cp -r /tmp/spack-staging/lib $SPACK_ROOT/lib                        \
 && cp -r /tmp/spack-staging/share $SPACK_ROOT/share                    \
 && cp -r /tmp/spack-staging/var $SPACK_ROOT/var                        \
 && cp -r /tmp/spack-staging/.git $SPACK_ROOT/.git                      \
 && rm -rf /tmp/spack-staging                                           \
 && echo 'export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH'\ 
        >> $SPACK_ROOT/share/setup-env.sh                               \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash                \
          /usr/sbin/docker-shell                                        \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash                \
          /usr/sbin/interactive-shell                                   \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash                \
          /usr/sbin/spack-env                                           \
 && echo "Part 2: Set target to generic x86_64"                         \
 && echo "packages:" > $SPACK_ROOT/etc/spack/packages.yaml              \
 && echo "  all:" >> $SPACK_ROOT/etc/spack/packages.yaml                \
 && echo "    target: [x86_64]" >> $SPACK_ROOT/etc/spack/packages.yaml  \
 && cat $SPACK_ROOT/etc/spack/packages.yaml                             \
 && echo "Part 3: Set config to allow use of more cores for builds"     \
 && echo "(and some other settings)"                                    \
 && echo "config:" > $SPACK_ROOT/etc/spack/config.yaml                  \
 && echo "  suppress_gpg_warnings: true"                                \
        >> $SPACK_ROOT/etc/spack/config.yaml                            \
 && echo "  build_jobs: 64" >> $SPACK_ROOT/etc/spack/config.yaml        \
 && echo "  install_tree:" >> $SPACK_ROOT/etc/spack/config.yaml         \
 && echo "    root: /opt/software" >> $SPACK_ROOT/etc/spack/config.yaml \
 && cat $SPACK_ROOT/etc/spack/config.yaml

SHELL ["docker-shell"]

## Setup spack buildcache mirrors, including an internal
## spack mirror using the docker build cache, and
## a backup mirror on the internal B010 network
RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    export PATH=$PATH:$SPACK_ROOT/bin                                   \
 && wget 10.10.241.24/spack-mirror/sodium.pub --no-check-certificate    \
 && spack gpg trust sodium.pub                                          \
 && spack mirror add silicon http://10.10.241.24/spack-mirror           \
 && spack mirror add docker /var/cache/spack-mirror                     \
 && spack mirror list

## Setup our custom environment and package overrides
COPY spack $SPACK_ROOT/eic-spack
RUN echo "repos:" > $SPACK_ROOT/etc/spack/repos.yaml                    \
 && echo " - $SPACK_ROOT/eic-spack" >> $SPACK_ROOT/etc/spack/repos.yaml \
 && mkdir /opt/spack-environment                                        \
 && mv $SPACK_ROOT/eic-spack/spack.yaml /opt/spack-environment/spack.yaml

## This variable will change whenevery either spack.yaml or our spack package
## overrides change, triggering a rebuild
ARG CACHE_BUST="hash"
## Now execute the main build (or fetch from cache if possible)
## note, no-check-signature is needed to allow the quicker signature-less
## packages from the internal (docker) buildcache
RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    cd /opt/spack-environment                                           \
 && ls /var/cache/spack-mirror                                          \
 && rm -r /usr/local                                                    \
 && spack env activate .                                                \
 && spack install -j64 --no-check-signature                             \
 && spack clean -a                                                  

## Update the local build cache if needed. Consists of 3 steps:
## 1. Remove the B010 network buildcache (silicon)
## 2. Get a list of all packages, and compare with what is already on
##    the buildcache (using package hash)
## 3. Add packages that need to be added to buildcache if any
RUN --mount=type=cache,target=/var/cache/spack-mirror                   \
    spack mirror remove silicon                                         \
 && spack buildcache list --allarch --long                              \
     | grep -v -e '---'                                                 \
     | sed "s/@.\+//"                                                   \
     | sort > tmp.buildcache.txt                                        \
 && spack find --no-groups --long                                       \
     | tail -n +2                                                       \
     | grep -v "==>"                                                    \
     | sed "s/@.\+//"                                                   \
     | sort > tmp.manifest.txt                                          \
 && comm -23 tmp.manifest.txt tmp.buildcache.txt                        \
     > tmp.needsupdating.txt                                            \
 && if [ $(wc -l < tmp.needsupdating.txt) -ge 1 ]; then                 \
     cat tmp.needsupdating.txt                                          \
        | awk '{print($2);}'                                            \
        | tr '\n' ' '                                                   \
        | xargs spack buildcache create -uaf -d /var/cache/spack-mirror \
     && spack buildcache update-index -d /var/cache/spack-mirror;       \
    fi                                                                  \
 && rm tmp.manifest.txt                                                 \
 && rm tmp.buildcache.txt                                               \
 && rm tmp.needsupdating.txt

## extra post-spack steps
## Including some small fixes:
RUN cd /opt/spack-environment                                           \
 && echo -n ""                                                          \
 && echo "Grabbing environment info"                                    \
 && spack env activate --sh -d .                                        \
        | sed "s?LD_LIBRARY_PATH=?&/lib/x86_64-linux-gnu:?"             \
        | sed '/MANPATH/ s/;$/:;/'                                      \
    > /etc/profile.d/z10_spack_environment.sh                           \
 && cd /opt/spack-environment && spack env activate .                   \
 && echo -n ""                                                          \
 && echo "Installing additional python packages"                        \
 && pip install --trusted-host pypi.org                                 \
                --trusted-host files.pythonhosted.org                   \
                --no-cache-dir                                          \
        ipython matplotlib scipy yapf pandas pycairo pyyaml lxml        \
        jupyter jupyterlab uproot pyunfold seaborn stashcp awkward

## make sure we have the entrypoints setup correctly
ENTRYPOINT []
CMD ["bash", "--rcfile", "/etc/profile", "-l"]
USER 0
WORKDIR /

## ========================================================================================
## STAGE 3: staging image with unnecessariy packages removed and stripped binaries
## ========================================================================================
FROM builder as staging

RUN cd /opt/spack-environment && spack env activate . && spack gc -y
# Strip all the binaries
# This reduces the image by factor of x2, so worth the effort
# note that we do not strip python libraries as can cause issues in some cases
RUN find -L /usr/local/*                                                \
         -type d -name site-packages -prune -false -o                   \
         -type f -not -name "zdll.lib"                                  \
         -exec realpath '{}' \;                                      \
      | xargs file -i                                                   \
      | grep 'charset=binary'                                           \
      | grep 'x-executable\|x-archive\|x-sharedlib'                     \
      | awk -F: '{print $1}' | xargs strip -s

## Address eic_container Issue #72
## missing precompiled headers for cppyy due to missing symlink in root
## install (should really be addressed by ROOT spack package)
RUN cd /opt/spack-environment && spack env activate .                   \
 && if [ ! -e $(spack location -i root)/lib/cppyy_backend/etc ]; then   \
      ln -sf $(spack location -i root)/etc                              \
             $(spack location -i root)/lib/cppyy_backend/etc;           \
    fi

RUN spack debug report                                                  \
      | sed "s/^/ - /" | sed "s/\* \*\*//" | sed "s/\*\*//"             \
    >> /etc/lager_info                                                  \
 && spack find --no-groups --long --variants | sed "s/^/ - /" >> /etc/lager_info

COPY lager-shell /usr/local/bin/lager-shell
COPY lager-info /usr/local/bin/lager-info
COPY entrypoint.sh /usr/local/sbin/entrypoint.sh
COPY lager-env.sh /etc/lager-env.sh
COPY profile.d/a00_cleanup.sh /etc/profile.d
COPY profile.d/z11_lager_env.sh /etc/profile.d
COPY singularity.d /.singularity.d

## Add minio client into /usr/local/bin
ADD https://dl.min.io/client/mc/release/linux-amd64/mc /usr/local/bin
RUN chmod a+x /usr/local/bin/mc

## random cleanup that's somehow necessary
RUN rm -rf /etc/profile.d/debuginfod.*

## ========================================================================================
## STAGE 4: Release image
## ========================================================================================
FROM base as release

LABEL maintainer="Sylvester Joosten <sjoosten@anl.gov>" \
      name="lager" \
      march="amd64"

## copy over everything we need from staging in a single layer :-)
RUN --mount=from=staging,target=/staging                                \
    rm -rf /usr/local                                                   \
 && cp -r /staging/opt/software /opt/software                           \
 && cp -r /staging/usr/._local /usr/._local                             \
 && cd /usr/._local                                                     \
 && PREFIX_PATH=$(realpath $(ls | tail -n1))                            \
 && echo "Found spack true prefix path to be $PREFIX_PATH"              \
 && cd -                                                                \
 && ln -s ${PREFIX_PATH} /usr/local                                     \
 && cp /staging/etc/profile.d/*.sh /etc/profile.d/                      \
 && cp /staging/etc/lager-env.sh /etc/lager-env.sh                      \
 && cp /staging/etc/lager_info /etc/lager_info                          \
 && cp -r /staging/.singularity.d /.singularity.d                        

## make sure we have the entrypoints setup correctly
ENTRYPOINT ["/usr/local/sbin/entrypoint.sh"]
CMD ["bash", "--rcfile", "/etc/profile", "-l"]
USER 0
WORKDIR /
SHELL ["/usr/local/bin/lager-shell"]

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
