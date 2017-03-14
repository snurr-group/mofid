.PHONY: all backup test diff ob_changes.patch

all:
	@echo "Sample make file for experimentation.  Still needs work.  Only backup implemented"

backup:
	rsync -av --exclude=".*" --exclude="openbabel/" --delete . ../../Box\ Sync/Projects/GitBackups/mofid

# Make this generic later on...
bin/sbu: src/sbu.cpp /usr/local/lib/openbabel/2.3.90/cifformat.so
	cd bin && cmake ../src && make

# Be careful: multi-line, nonescaped commands in Make run in separate shells
# Generic rules for compiling relevant (modified by me) formats
/usr/local/lib/openbabel/2.3.90/cifformat.so: openbabel/src/formats/cifformat.cpp openbabel/src/formats/systreformat.cpp
	cd /cygdrive/c/Users/Benjamin/Git/mofid/openbabel/build; \
	make cifformat; \
	make systreformat; \
	make install/fast
	dos2unix ${BABEL_DATADIR}/*.txt  # Necessary so C++ can actually read the files

diff: ob_changes.patch

ob_changes.patch:
	git diff --no-prefix 7810ca7bb1beef14b2a62cf5bad3a8551b187824 -- openbabel/*.cpp openbabel/*.h ':!openbabel/data/*' > $@
	# Lists my changes to the main OpenBabel code
