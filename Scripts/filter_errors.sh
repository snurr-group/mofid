#!/bin/bash
# Filters output json results to focus on errors.
# Only shows problematic MOFs from the JSON of check_mof_linkers.py
# instead of manually sorting through the output file.

JQ_EXE="Resources/External/jq"
STDIN_CONTENTS=$(cat)

echo "$STDIN_CONTENTS" | "$JQ_EXE" '.mofs | map(select((.errors | length) > 0))'


# Other information:
# jq is inspired by unix pipes for commands and args.
# See documentation at https://stedolan.github.io/jq/manual/
# Another example query to filter out the longest simulations:
# cat minimal.json | Resources/External/jq '.mofs | map(select(.time >= 10)) | sort_by(.time) | reverse'

# Alternatives considered before implementing jq:
#* [jsonpath](https://goessner.net/articles/JsonPath/) looks like a great standard that does everything I want
#	* Unfortunately, the [Python implementation](https://github.com/kennknowles/python-jsonpath-rw) doesn't appear to include selector queries
#	* And the -ext version doesn't look viable, either
#	* Maybe there would be another useful implementation: <https://www.pluralsight.com/blog/tutorials/introduction-to-jsonpath>
#* [ObjectPath](http://objectpath.org/) appears deprecated but otherwise exactly what I want
#* Apparently XQuery3.1 will have native JSON support, so hopefully this will not be a problem much longer.

