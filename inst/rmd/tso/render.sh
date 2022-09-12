#!/usr/bin/env bash

here="/Users/pdiakumis/projects/dracarys"
dd="${here}/nogit/tso"
pref="SBJ02096/PRJ221878_L2201128"
prefix="${dd}/${pref}"

quarto render \
    "tso.qmd" \
    -P prefix:${prefix} \
    --to html \
    --output "SBJ02096_dracarys.html"
