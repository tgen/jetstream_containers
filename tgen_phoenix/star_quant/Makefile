.PHONY: build publish clean

MAKEFILE_PATH           := $(abspath $(lastword $(MAKEFILE_LIST)))
IMAGE_NAME              ?= tgen_phoenix_$$(basename $$(dirname ${MAKEFILE_PATH}) )
SOURCE_CONTAINER_REGISTRY  ?= localhost:5000
PUBLISH_CONTAINER_REGISTRY ?= localhost:5000
IMAGE_PIPELINE_VERSION  ?= 20220503
IMAGE_ARCH              ?= skylake
IMAGE_TAG               ?= ${PUBLISH_CONTAINER_REGISTRY}/${IMAGE_NAME}:${IMAGE_PIPELINE_VERSION}-${IMAGE_ARCH}
LOCAL_SOURCES_SERVER    ?= localhost:8000

default: build publish

build: Dockerfile
	docker build \
	  --network=host \
	  --build-arg SOURCE_CONTAINER_REGISTRY=${SOURCE_CONTAINER_REGISTRY} \
	  --build-arg IMAGE_ARCH=${IMAGE_ARCH} \
	  --build-arg LOCAL_SOURCES_SERVER=${LOCAL_SOURCES_SERVER} \
	  -t ${IMAGE_TAG}  \
	  .
	  
publish:
	docker push ${IMAGE_TAG} 

clean:
	docker build rm ${IMAGE_TAG} 
