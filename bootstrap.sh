#!/bin/sh

# ensure m4 dir exists 
mkdir -p m4
# ensure build-aux exists
mkdir -p build-aux

autoreconf -s -i -m -f
