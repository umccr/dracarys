date_start="2024-08-29"
date_end="2024-09-07"
out="sash_${date_start}_${date_end}.html"

quarto render summary_sash.qmd \
    -P date_start:${date_start} \
    -P date_end:${date_end} \
    -o ${out} \
    --output-dir nogit/html
