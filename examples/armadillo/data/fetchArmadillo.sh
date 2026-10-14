#!/bin/bash
#
# Make sure you chmod +x fetchArmadillo.sh to make it executable 
#
# Fetch the Stanford Armadillo geometry from the github repository
URL="https://github.com/alecjacobson/common-3d-test-models/raw/master/data/armadillo.obj"
# Download the geometry
curl -L -o armadillo.obj "$URL"
