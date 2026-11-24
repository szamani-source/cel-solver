# celFoam in a container: works on any Linux with Docker or Podman (INSTALL.md, route C).
#
#   docker build -t celfoam:1.1.0 .                 # or: ./install.sh --container
#   docker run --rm -it -v "$PWD":/work celfoam:1.1.0            # shell, environment loaded
#   docker run --rm -v "$PWD":/work celfoam:1.1.0 ./Allrun       # one command in /work
#
# Base: the OpenFOAM v2512 image that CI and tools/setup_openfoam.sh use, pinned to the
# same amd64 manifest (decision 0006), so the container has the same OpenFOAM bits.
FROM opencfd/openfoam-dev:2512@sha256:0c279c90cc0531e4f5acfe5365c529d5794a50fc7d8e4a0a45037c5f404f21b9

USER root
# The image's apt source for dl.openfoam.com is not needed (OpenFOAM is already in the
# image) and has failed during mirror syncs, so it is dropped before the update (as CI).
RUN grep -rl openfoam /etc/apt/sources.list.d/ 2>/dev/null | xargs -r rm -f \
 && apt-get update -qq \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq --no-install-recommends \
        python3 python3-numpy python3-matplotlib \
 && rm -rf /var/lib/apt/lists/*

COPY . /tmp/celfoam-src
RUN /tmp/celfoam-src/install.sh --prefix /opt/celfoam \
        --openfoam /usr/lib/openfoam/openfoam2512 --no-link --no-test \
 && ln -s /opt/celfoam/bin/celfoam /usr/local/bin/celfoam \
 && rm -rf /tmp/celfoam-src \
 && chmod -R a+rX /opt/celfoam

ENV CELFOAM_DIR=/opt/celfoam OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /work
ENTRYPOINT ["/opt/celfoam/bin/celfoam"]
