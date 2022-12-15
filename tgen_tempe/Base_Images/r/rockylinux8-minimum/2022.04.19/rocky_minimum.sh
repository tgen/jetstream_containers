#!/bin/bash

# using date to set a version for.
# To maintain different base images with different build dates
# Also for informational purposes.
buildDate=$(date '+%Y.%m.%d')

# Create an empty image with [scratch]
newcontainer=$(buildah from scratch)
# Mount the scratch container
scratchmnt=$(buildah mount $newcontainer)
# install packages to [scratch] container
dnf -y group install "Minimal Install" --releasever 8 --installroot=$scratchmnt

# unmount
buildah umount $newcontainer

# add image
buildah commit $newcontainer rockylinux8-minimum:${buildDate}
