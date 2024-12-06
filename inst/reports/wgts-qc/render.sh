date1="${1:-$(date +%Y-%m-%d)}"
out="alignqc_${date1}.html"
tidy_data="nogit/tidy_data_rds/${date1}_wgts.rds"

quarto render summary.qmd \
    -P tidy_data:${tidy_data} \
    -o ${out} \
    --output-dir "nogit/html"
