#!/bin/sh
if [ ! -d src/liblsoda ]; then
    git clone https://github.com/richfitz/lsoda src/liblsoda
else
    git -C src/liblsoda pull
fi
