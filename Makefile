.PHONY: all backup test diff ob_changes.patch

all:
	@echo "Sample make file for experimentation.  Still needs work.  Only backup implemented"

backup:
	rsync -av --exclude=".*" --exclude="openbabel/" --delete . ../../Box\ Sync/Projects/GitBackups/mofid

# Make this generic later on...
bin/sbu: src/sbu.cpp openbabel/build/lib/cifformat.so
	cd bin && make

# Be careful: multi-line, nonescaped commands in Make run in separate shells
# Generic rules for compiling relevant (modified by me) formats
openbabel/build/lib/cifformat.so: openbabel/src/formats/cifformat.cpp openbabel/src/formats/systreformat.cpp
	cd mofid/openbabel/build; \
	make cifformat; \
	make systreformat; \
	make install/fast

diff: ob_changes.patch

ob_changes.patch:
	git diff --no-prefix 7810ca7bb1beef14b2a62cf5bad3a8551b187824 -- openbabel/*.cpp openbabel/*.h ':!openbabel/data/*' > $@
	# Lists my changes to the main OpenBabel code

setup:
	cd openbabel; \
	mkdir build installed; \
	cd build; \
	cmake -DCMAKE_INSTALL_PREFIX=../installed ..; \
	make -j2;  # Build in parallel \
	make install; # This is apparently a crucial step to get the sbu linker working correctly \
	cd ../../; \
	mkdir bin; \
	cd bin; \
	cmake -DOpenBabel2_DIR=../openbabel/build ../src/; \
	make
	# Sets up all the cmake details, so that usage is as simple as
	# `bin/sbu MOF.cif` and re-compilation is as easy as `make bin/sbu`
