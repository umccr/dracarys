#!/usr/bin/env bash

here="/Users/pdiakumis/projects/dracarys"
dd="${here}/nogit/tso/2022-12-13"
pref="SBJ00757/dracarys_gds_sync/PRJ222072_L2201397"
prefix="${dd}/${pref}"

quarto render \
    "tso.qmd" \
    -P prefix:${prefix} \
    --to html \
    --output "SBJ00757_dracarys.html"
