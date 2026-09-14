#!/usr/bin/env sh
# Assemble the static stlite site into dist/ (used by CI and for local preview:
#   sh scripts/build_site.sh && python -m http.server -d dist 8000)
set -eu
cd "$(dirname "$0")/.."
rm -rf dist
mkdir -p dist/transquant_w
cp web/index.html app.py dist/
cp src/transquant_w/__init__.py src/transquant_w/core.py src/transquant_w/ui.py dist/transquant_w/
cp -r static dist/
touch dist/.nojekyll
echo "built dist/"
