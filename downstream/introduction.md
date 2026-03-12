## Input

Pipeline nay can file variant-level theo mau:

    variant_concordance_10x.tsv
    variant_concordance_20x.tsv
    variant_concordance_30x.tsv
    variant_concordance_50x.tsv

Va moi file phai co cot:

    variant_id    truth    HC    DV    ST    FB    DS


## Scripts

```bash
# pwd: variant-calling-benchmark
export CONCORDANCE_SEARCH_ROOT="results/eval"
export CONCORDANCE_DIR="results/eval/concordance"

python downstream/build_concordance_from_vcfs.py \
    --project-root . \
    --search-root "$CONCORDANCE_SEARCH_ROOT" \
    --out-dir "$CONCORDANCE_DIR"
python downstream/extract_error_patterns.py --project-root . --vcfeval-dir "$CONCORDANCE_DIR"
python downstream/annotate_variants_strata.py --project-root .
```

## Output

- results/analysis/error_patterns
- results/analysis/stratification
- results/analysis/summary_tables
- results/plots/downstream
