#!/usr/bin/env bash

set -euo pipefail

gatk_ver=4.2.1.0
gatk_url=https://github.com/broadinstitute/gatk/releases/download/${gatk_ver}/gatk-${gatk_ver}.zip

wget ${gatk_url}
unzip -q gatk-${gatk_ver}.zip
rm gatk-${gatk_ver}.zip
gatk-${gatk_ver}/gatk --version

>&2 echo Done
exit 0

