#!/bin/bash

# kinda like `quarto render`, but in order

OUTPUT_DIR="../docs/"
TARGETS=(
  index.md
  installation.md

  tutorial-fmri-spm/index.md
  tutorial-fmri-spm/preparation.qmd
  tutorial-fmri-spm/model.qmd
  tutorial-fmri-spm/cvmanova.qmd
  tutorial-fmri-spm/pairwise.qmd
  tutorial-fmri-spm/cvcrossmanova.qmd

  tutorial-other/index.md
  tutorial-other/preparation.qmd

  reference.qmd
)


shopt -s dotglob    # include hidden files in `*`
shopt -s globstar   # enable `**`
shopt -s nullglob   # no matches expands to nothing
set -e              # abort immediately upon SIGINT


echo
echo "––– Update bibliography from Zotero ––––––––––––––––––––––––––––––––––––"
echo
curl -s "http://127.0.0.1:23119/better-bibtex/export/collection?/1/cvcrossmanova.yaml" -o bibliography.yaml
if [ $? -ne 0 ]; then
  echo "Warning: Zotero does not appear to be running."
fi

echo
echo "––– Deleting caches and output -––––––––––––––––––––––––––––––––––––––––"
echo

names=("${OUTPUT_DIR%/}"/* **/.quarto/ **/.jupyter_cache/)
printf "  %s\n" "${names[@]}"
read -p "delete? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
  rm -r "${names[@]}"
fi


echo
echo "––– Rendering with Quarto ––––––––––––––––––––––––––––––––––––––––––––––"
echo

# "quarto render" command. This version runs the command in a micromamba environment "std" in which the the MKernel Jupyter kernel is installed.
Q="micromamba run -n std quarto render"

for T in ${TARGETS[@]}
do
  $Q "$T"
done

echo
echo "––– Transforming text to paths in SVGs –––––––––––––––––––––––––––––––––"
echo

for F in "${OUTPUT_DIR%/}"/**/*.svg
do
  echo "processing '$F'"
  inkscape --export-type=svg --export-overwrite --export-text-to-path "$F" \
    |& grep -E -v "WARNING|appmenu-gtk-module|^$" \
    || true
  echo
done


echo
