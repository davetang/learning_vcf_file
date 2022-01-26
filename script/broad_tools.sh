#!/usr/bin/env bash
#
# Download GATK, Cromwell, and WOMtool
#

set -euo pipefail

# GATK
gatk_ver=4.2.0.0
gatk_url=https://github.com/broadinstitute/gatk/releases/download/${gatk_ver}/gatk-${gatk_ver}.zip

if [[ ! -d gatk-${gatk_ver} ]]; then
   wget ${gatk_url}
   unzip -q gatk-${gatk_ver}.zip
   rm gatk-${gatk_ver}.zip
   gatk-${gatk_ver}/gatk --version
else
   >&2 echo gatk-${gatk_ver} already exists
fi

# Cromwell
cromwell_ver=65
cromwell_url=https://github.com/broadinstitute/cromwell/releases/download/${cromwell_ver}/cromwell-${cromwell_ver}.jar
if [[ ! -e cromwell-${cromwell_ver}.jar ]]; then
   wget ${cromwell_url}
   java -jar cromwell-${cromwell_ver}.jar --version
else
   >&2 echo cromwell-${cromwell_ver}.jar already exists
fi

# WOMtool
womtool_ver=65
womtool_url=https://github.com/broadinstitute/cromwell/releases/download/${womtool_ver}/womtool-${womtool_ver}.jar
if [[ ! -e womtool-${womtool_ver}.jar ]]; then
   wget ${womtool_url}
   java -jar womtool-${womtool_ver}.jar --version
else
   >&2 echo womtool-${womtool_ver}.jar already exists
fi

