ARG SOURCE_DOCKER_REGISTRY=localhost:5000
ARG IMAGE_ARCH=sandybridge

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_python-3.8.6.headless:${IMAGE_ARCH} as opt_python

FROM ${SOURCE_DOCKER_REGISTRY}/leanix-1.0_opt_base:${IMAGE_ARCH}

COPY --from=opt_python /opt/ /opt/
COPY capture_kits.csv /opt/phoenix/
COPY make_exome_refpack.py /opt/utility_scripts/
COPY formats /opt/lib/python3.8/site-packages/formats

LABEL org.opencontainers.image.source=https://github.com/tgen/jetstream_containers

