# Shell script to create a minimum debian 10 image that is used to build debian10base.
# Use docker file to add the LABEL to this image.

# Create an empty image with [scratch]
newcontainer=$(buildah from scratch)
# Moun the scratch container
scratchmnt=$(buildah mount $newcontainer)
# Install base debian image with debootstrap
dnf install debootstrap
mount -i -o remount,exec,dev $scratchmnt
debootstrap --variant=minbase --arch amd64 buster $scratchmnt http://http.debian.net/debian

# unmount
buildah umount $newcontainer
# Add image
buildah commit $newcontainer debian10-minimum:latest
