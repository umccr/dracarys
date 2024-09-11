date_start="2024-08-29"
date_end="2024-09-01"
out="umccrise_${date_start}_${date_end}.html"

quarto render summary_umccrise.qmd \
    -P date_start:${date_start} \
    -P date_end:${date_end} \
    -o ${out} \
    --output-dir nogit/html
