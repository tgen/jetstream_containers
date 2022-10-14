# Tgen-Pipelines

Open Container Initiative containers for use at TGen.

## Documentation Pending
---
## Objective: 
- Create containers with only providing a Dockerfile

### Methods:

    Will use folder path to name and tag containers
    Example:

    python/2.7/Dockerfile => <URL>/python:2.7

### 1.) Create a New Container 

- git pull <url> 
- create a folder/use existing folder with the image name 
- create a folder and name it with the docker tag
- create dockerfile 
- git add . / git commit -m "" / git push
- YOUR CONTAINER HAS BEEN BUILT
- once container is working properly submit a MR to add it to the master branch and protect it

### 2.) Update an Existing Container 

- git pull <url> 
- use existing folder with the image name/tag
- update dockerfile 
- git add . / git commit -m "" / git push
- submit MR from dev to master
- await for MR approval
- Once MR has been approved, it will trigger the pipeline and overwrite the existing dockerfile with the new configurations

---
## How to create an image from a Dockerfile

    1.) sudo podman login --tls-verify=false pbc-art-prd01.ad.tgen.org

    2.) sudo buildah bud -f {dockerfile)} -t {image name} .

    3.) sudo podman tag {image:tag} pbc-art-prd01.ad.tgen.org/hpc-local-containers/{image name}:{tag}

    4.) sudo podman push pbc-art-prd01.ad.tgen.org/hpc-local-containers/{image name}:{tag}

## How to pull an image from Artifactory 
    1.) sudo podman login --tls-verify=false pbc-art-prd01.ad.tgen.org

    2.) sudo podman pull pbc-art-prd01.ad.tgen.org/hpc-virtual-containers/{image name}:{tag}
## 