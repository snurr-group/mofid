#!/bin/bash
# Downloads jq program.  Used internally by base Makefile

JQ_BIN="UNKNOWN"
kernel=$(uname -s)

if [[ $kernel == *NT* ]]
then
	JQ_BIN="jq-win64.exe"
elif [[ $kernel == *Linux* || $kernel == *linux* ]]
then
	JQ_BIN="jq-linux64"
elif [[ $kernel == *Darwin* ]]
then
	JQ_BIN="jq-osx-amd64"
else
	echo "Unknown platform.  Detected: $kernel"
	exit 3
fi

wget "https://github.com/stedolan/jq/releases/download/jq-1.6/${JQ_BIN}"
wget "https://raw.githubusercontent.com/stedolan/jq/master/sig/v1.6/sha256sum.txt"

mv sha256sum.txt check_jq_sha256.txt
echo "Verifying checksum integrity of jq executable:"
grep $JQ_BIN check_jq_sha256.txt | sha256sum --check || exit 55
rm check_jq_sha256.txt

mv $JQ_BIN jq
chmod +x jq

