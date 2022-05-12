.PHONY: build push

MAKEFILE_PATH           := $(abspath $(lastword $(MAKEFILE_LIST)))
IMAGE_NAME              ?= tgen_phoenix_$$(basename $$(dirname ${MAKEFILE_PATH}) )
SOURCE_CONTAINER_REGISTRY  ?= localhost:5000
PUBLISH_CONTAINER_REGISTRY ?= localhost:5000
IMAGE_PIPELINE_VERSION  ?= 20220503
IMAGE_ARCH              ?= skylake
IMAGE_TAG               ?= ${PUBLISH_CONTAINER_REGISTRY}/${IMAGE_NAME}:${IMAGE_PIPELINE_VERSION}-${IMAGE_ARCH}

default: build push

build: Dockerfile
	buildah bud \
	    -t ${IMAGE_TAG} \
		-f ./Dockerfile \
		--layers=true \
		--pull \
		--network=host \
		--build-arg SOURCE_CONTAINER_REGISTRY=${SOURCE_CONTAINER_REGISTRY} \
		--build-arg IMAGE_ARCH=${IMAGE_ARCH}

push:
	buildah push ${IMAGE_TAG}
