#!/bin/bash

echo "render.sh $1"

shopt -s dotglob    # include hidden files in `*`
shopt -s globstar   # enable `**`
shopt -s nullglob   # no matches expand to nothing

# pre- and post-render
case $1 in

  pre)
    if [[ -n "$QUARTO_PROJECT_RENDER_ALL" ]]  # full render
    then
      echo "Deleting caches and output"
      names=("$QUARTO_PROJECT_OUTPUT_DIR"* **/.quarto/ **/.jupyter_cache/)
      printf "  %s\n" "${names[@]}"
      read -p "delete? " -n 1 -r
      echo
      if [[ $REPLY =~ ^[Yy]$ ]]
      then
          rm -r "${names[@]}"
      fi
    fi
    ;;

  post)
    if [[ -n "$QUARTO_PROJECT_RENDER_ALL" ]]  # full render
    then
      echo "Transforming text to paths in SVGs"
      for fn in "$QUARTO_PROJECT_OUTPUT_DIR"**/*.svg
      do
        echo "processing '$fn'"
        inkscape --export-type=svg --export-overwrite \
          --export-text-to-path "$fn" \
          |& grep -E -v "WARNING|appmenu-gtk-module|^$"
      done
    fi
    ;;

  *)
    echo -n "Only to be used from '_quarto.yml'!"
    ;;

esac
