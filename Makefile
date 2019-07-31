.PHONY: all backup test diff ob_changes.patch init eclipse web init-web github-web html one exe btc

all:
	@echo "Sample make file for experimentation.  Still needs work.  Only backup implemented"

backup:
	rsync -av --exclude=".*" --exclude="openbabel/" --delete . ../../Box\ Sync/Projects/GitBackups/mofid

# Make this generic later on...
bin/sbu: src/sbu.cpp openbabel/build/lib/cifformat.so
	cd bin && make sbu
bin/sobgrep: src/sobgrep.cpp openbabel/build/lib/cifformat.so
	cd bin && make sobgrep
bin/searchdb: src/searchdb.cpp openbabel/build/lib/cifformat.so
	cd bin && make searchdb
bin/tsfm_smiles: src/tsfm_smiles.cpp openbabel/build/lib/cifformat.so
	cd bin && make tsfm_smiles

exe:
	cd bin && make

one:
	cd bin && make; \
	cd ..; \
	bin/sbu Resources/TestCIFs/P1-IRMOF-1.cif

btc:
	bin/sbu Resources/TestCIFs/P1-Cu-BTC.cif

# Be careful: multi-line, nonescaped commands in Make run in separate shells
# Generic rules for compiling relevant (modified by me) formats
openbabel/build/lib/cifformat.so: openbabel/src/formats/cifformat.cpp openbabel/src/mol.cpp
	cd openbabel/build; \
	make cifformat; \
	make install/fast

diff: ob_changes.patch

ob_changes.patch:
	git diff --no-prefix 7810ca7bb1beef14b2a62cf5bad3a8551b187824 -- openbabel/*.cpp openbabel/*.h ':!openbabel/data/*' ':!openbabel/test/*' > $@
	# Lists my changes to the main OpenBabel code

test: bin/sbu
	python Python/check_mof_composition.py



init:
	cd openbabel; \
	mkdir build installed; \
	cd build; \
	cmake -DCMAKE_INSTALL_PREFIX=../installed -DBUILD_GUI=OFF -DEIGEN3_INCLUDE_DIR=../eigen ..; \
	make -j2 || exit 2; \
	make install; \
	cd ../../; \
	mkdir bin; \
	cd bin; \
	cmake -DOpenBabel2_DIR=../openbabel/build ../src/; \
	make
	# Sets up all the cmake details, so that usage is as simple as
	# `bin/sbu MOF.cif` and re-compilation is as easy as `make exe`

eclipse:
	cd bin; \
	cmake -G "Eclipse CDT4 - Unix Makefiles" ../src; \

# Emscripten web content below
# In my current Windows setup, these must all be run within Git Bash.
# Emscripten should be installed in /mofid/../emsdk/ using the instructions at
# https://emscripten.org/docs/getting_started/downloads.html
# Node.js must also be installed to successfully compile webGavrog.
#
# After having some problems in Windows, I switched to the quest Linux computing cluster.
# The steps are installing emscripten using the link above.  Then you need to load modules and run the make targets:
# module load clang/3.4.1 && module load gcc/6.4.0 && module load cmake/3.12.0
# make init-web && make web && make github-web

init-web:
	../emsdk/emsdk activate latest && source ../emsdk/emsdk_env.sh; \
	cd openbabel; \
	mkdir embuild eminstall; \
	cd embuild; \
	emcmake cmake .. -DCMAKE_INSTALL_PREFIX=../eminstall/ -DBUILD_GUI=OFF -DEIGEN3_INCLUDE_DIR=../eigen -DENABLE_TESTS=OFF -DBUILD_SHARED=OFF -DWITH_STATIC_INCHI=ON -DWITH_STATIC_XML=ON -DCMAKE_CXX_FLAGS="-O3 -s WASM=1"; \
	cd ../..; \
	mkdir embin; \
	cd embin; \
	emcmake cmake -DOpenBabel2_DIR=../openbabel/embuild -static ../src/ -DCMAKE_CXX_FLAGS="-O3 --preload-file ../openbabel/data@/ob_datadir/ --preload-file ../src/Web/web_data@/web_data/ --preload-file ../Resources/RCSRnets.arc@/RCSRnets.arc --pre-js ../src/pre_emscripten.js -s TOTAL_MEMORY=128MB -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS=\"['ccall', 'cwrap', 'UTF8ToString']\""; \
	mkdir kekule; \
	cd kekule; \
	unzip ../../Resources/kekule.release.0.7.5.170624.zip; \
	cd ..; \
	unzip ../Resources/webGavrog-20190721.zip && mv webGavrog-5649192c38329d8616ce8e2204787bdc945f6dc4 webGavrog-build; \
	mkdir webGavrog

openbabel/embuild/obabel.js:
	../emsdk/emsdk activate latest && source ../emsdk/emsdk_env.sh; \
	cd openbabel/embuild; \
	emmake make; \
	emmake make install

web: embin/sbu.js embin/webGavrog/main.js html

github-web: web
	rm -r ../web-mofid/*; \
	rm -f embin/*.wast; \
	cp embin/*.* embin/.gitignore ../web-mofid; \
	cp -r embin/webGavrog embin/kekule ../web-mofid;
	# Commit and push the changes in the web-mofid repo to update
	# the live website on Github.

html: src/Web/*.* src/Web/.gitignore Resources/ngl.js
	cp $^ embin/

embin/webGavrog/main.js: src/Web/gavrog_override/*.js
	cp -r src/Web/gavrog_override/* embin/webGavrog-build/; \
	cd embin/webGavrog-build; \
	npm install; \
	npm run build; \
	cp dist/*.js ../webGavrog/
	# If better debug info is needed for webGavrog, replace `npm run build` with `npm run build-dev`

embin/sbu.js: src/sbu.cpp openbabel/embuild/obabel.js src/pre_emscripten.js
	../emsdk/emsdk activate latest && source ../emsdk/emsdk_env.sh; \
	cd embin; \
	emmake make

