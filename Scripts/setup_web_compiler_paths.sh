#!/bin/bash
# Sets up paths for emscripten in the current shell.
# Requires emscripten to be installed.
# Run this program via source import_emscripten.sh
# TODO: Consider using a variable like $HOME for the paths for portability

export PATH=/c/Users/Benjamin/Programs/emsdk-portable:/c/Users/Benjamin/Programs/emsdk-portable/clang/e1.37.18_64bit:/c/Users/Benjamin/Programs/emsdk-portable/node/4.1.1_64bit/bin:/c/Users/Benjamin/Programs/emsdk-portable/python/2.7.5.3_64bit:/c/Users/Benjamin/Programs/emsdk-portable/java/7.45_64bit/bin:/c/Users/Benjamin/Programs/emsdk-portable/emscripten/1.37.18:$PATH
export EMSDK=c:/Users/Benjamin/Programs/emsdk-portable
export EM_CONFIG=/c/Users/Benjamin/.emscripten
export BINARYEN_ROOT=/c/Users/Benjamin/Programs/emsdk-portable/clang/e1.37.18_64bit/binaryen
export JAVA_HOME=/c/Users/Benjamin/Programs/emsdk-portable/java/7.45_64bit
export EMSCRIPTEN=/c/Users/Benjamin/Programs/emsdk-portable/emscripten/1.37.18
emsdk activate latest
