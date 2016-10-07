#!/bin/sh
if [ ! -d src/liblsoda ]; then
    git clone https://github.com/richfitz/liblsoda src/liblsoda
else
    git -C src/liblsoda pull
fi
