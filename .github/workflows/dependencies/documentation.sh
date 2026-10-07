#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Andrew Myers

set -eu -o pipefail

sudo apt-get update

# The docs job only builds HTML with Sphinx; doxygen and LaTeX are not used.
sudo apt-get install -y --no-install-recommends\
    build-essential \
    pandoc

