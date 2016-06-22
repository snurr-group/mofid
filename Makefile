.PHONY: all backup test

all:
	@echo "Sample make file for experimentation.  Still needs work.  Only backup implemented"

backup:
	rsync -av --exclude=".*" --delete . ../../Box\ Sync/Projects/GitBackups/mofid

