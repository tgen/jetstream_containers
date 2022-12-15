# dnf update or another command to update repositories
# needs to be run on images created with this base image
# to ensure that the epel and powertools repos are added.

FROM ghcr.io/tgen/jetstream_containers/rockylinux8-minimum:2022.04.19

RUN dnf install -y rsync dnf-plugins-core epel-release procps
RUN dnf config-manager --set-enabled powertools

# ghcr.io/tgen/jetstream_containers/rocky8_with_epel:2022.04.19