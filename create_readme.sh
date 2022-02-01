#!/usr/bin/env bash

set -euo pipefail

out_md=tmp.md
Rscript -e "rmarkdown::render('readme.Rmd', output_file=\"$out_md\")"
cp -f $out_md mkdocs/docs/index.md
tools/bin/gh-md-toc $out_md > toc

cat toc <(echo) <(date) <(echo) $out_md > README.md

>&2 echo Done!
exit 0

