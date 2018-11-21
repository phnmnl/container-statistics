#!/bin/bash

# install required packages
apt-get update -y && apt-get install -y --no-install-recommends wget ca-certificates

# retrieve test case data
mkdir /tmp/testfiles
wget -O /tmp/testfiles/peaktable.txt https://github.com/phnmnl/container-statistics/raw/develop/testfiles/peaktable.txt
wget -O /tmp/testfiles/var.txt https://github.com/phnmnl/container-statistics/raw/develop/testfiles/var.txt
wget -O /tmp/testfiles/metadata.txt https://github.com/phnmnl/container-statistics/raw/develop/testfiles/metadata.txt

# perform test case
/usr/local/bin/univeriateLimma.r -in input=/tmp/testfiles/peaktable.txt -s /tmp/testfiles/metadata.txt -p /tmp/testfiles/var.txt -out /tmp/test.txt -g Groups --contrasts "A-B"

# check whether files were created
if [ ! -f /tmp/test.txt ]; then exit 1; fi

# remove files
rm /tmp/test.txt
