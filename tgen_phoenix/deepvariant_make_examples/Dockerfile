# ARG SOURCE_CONTAINER_REGISTRY=localhost:5000
# ARG IMAGE_ARCH=skylake

# FROM ${SOURCE_CONTAINER_REGISTRY}/ubuntu-16.04_opt_deepvariant-0.10.0:${IMAGE_ARCH} 

# ENV PATH /opt/deepvariant/bin/:/opt/bin/:${PATH}
# ENV LD_LIBRARY_PATH /opt/lib64/:/opt/lib/:${LD_LIBRARY_PATH}

# RUN mkdir -p /usr/bin/ && echo "#!/bin/bash" > /usr/bin/module && chmod a+x /usr/bin/module

FROM ghcr.io/tgen/jetstream_containers/tgen_phoenix_deepvariant_make_examples:20210113-skylake