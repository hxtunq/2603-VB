# Downstream Run Quick (Ubuntu Bash)

## 0) Luu y ve input concordance

Pipeline nay can file variant-level theo mau:

    variant_concordance_10x.tsv
    variant_concordance_20x.tsv
    variant_concordance_30x.tsv
    variant_concordance_50x.tsv

Va moi file phai co cot:

    variant_id    truth    HC    DV    ST    FB    DS

Neu ban chi co `concordance_matrix.tsv` dang pairwise summary
(`Coverage, CallerA, CallerB, Concordance`) thi KHONG du thong tin
de trich all-5-FN/all-5-FP/single-caller-FNFP o muc bien the.

## 1) Chay tung buoc (khong AlphaGenome)

    cd /home/hxt/phase-1/2603-VB

        # Neu CHUA co variant_concordance_<cov>x.tsv, script se tu build tu tp/fp/fn
        export CONCORDANCE_SEARCH_ROOT="results/eval"
        export CONCORDANCE_DIR="results/eval/concordance"

        python downstream/build_concordance_from_vcfs.py \
            --project-root . \
            --search-root "$CONCORDANCE_SEARCH_ROOT" \
            --out-dir "$CONCORDANCE_DIR"

    python downstream/extract_error_patterns.py --project-root . --vcfeval-dir "$CONCORDANCE_DIR"
    python downstream/annotate_variants_strata.py --project-root .
    python downstream/aggregate_alphagenome_scores.py --project-root .

    Rscript downstream/upset_plot_enhanced.R .
    Rscript downstream/plot_strata_characterization.R .
    Rscript downstream/analyze_coverage_sensitivity.R . "./$CONCORDANCE_DIR"

## 2) Chay AlphaGenome roi merge

    cd /home/hxt/phase-1/2603-VB
    export ALPHAGENOME_API_KEY="YOUR_KEY"

    bash downstream/run_alphagenome_batch.sh
    python downstream/aggregate_alphagenome_scores.py --project-root .

## 3) All-in-one (khong AlphaGenome)

    cd /home/hxt/phase-1/2603-VB
    export CONCORDANCE_SEARCH_ROOT="results/eval"
    export CONCORDANCE_DIR="results/eval/concordance"
    export AUTO_BUILD_CONCORDANCE=1
    bash downstream/run_downstream_pipeline.sh

## 4) All-in-one (co AlphaGenome)

    cd /home/hxt/phase-1/2603-VB
    export CONCORDANCE_SEARCH_ROOT="results/eval"
    export CONCORDANCE_DIR="results/eval/concordance"
    export AUTO_BUILD_CONCORDANCE=1
    export ALPHAGENOME_API_KEY="YOUR_KEY"
    export RUN_ALPHAGENOME=1
    bash downstream/run_downstream_pipeline.sh

## Output chinh

- results/analysis/error_patterns
- results/analysis/stratification
- results/analysis/alphagenome
- results/analysis/summary_tables
- results/plots/downstream
