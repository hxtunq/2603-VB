# Downstream Run Quick

## 1) Chay nhanh tung buoc (PowerShell)

    cd D:\Git-Page\variant-calling-benchmark

    python downstream\extract_error_patterns.py --project-root .
    python downstream\annotate_variants_strata.py --project-root .

    # Neu chua chay AlphaGenome thi van co the chay buoc merge (se ra canh bao nhe)
    python downstream\aggregate_alphagenome_scores.py --project-root .

    Rscript downstream\upset_plot_enhanced.R .
    Rscript downstream\plot_strata_characterization.R .
    Rscript downstream\analyze_coverage_sensitivity.R .

## 2) Neu muon chay AlphaGenome truoc khi merge

    cd D:\Git-Page\variant-calling-benchmark

    # Dat API key
    $env:ALPHAGENOME_API_KEY = "YOUR_KEY"

    # Chay batch score cho 5 callers x 4 coverages x FN/FP
    bash downstream/run_alphagenome_batch.sh

    # Merge lai diem vao bang downstream
    python downstream\aggregate_alphagenome_scores.py --project-root .

## 3) Lenh all-in-one (khong chay AlphaGenome)

    cd D:\Git-Page\variant-calling-benchmark
    bash downstream/run_downstream_pipeline.sh

## 4) Lenh all-in-one (co chay AlphaGenome)

    cd D:\Git-Page\variant-calling-benchmark
    $env:ALPHAGENOME_API_KEY = "YOUR_KEY"
    $env:RUN_ALPHAGENOME = "1"
    bash downstream/run_downstream_pipeline.sh

## Output chinh

- results/analysis/error_patterns
- results/analysis/stratification
- results/analysis/alphagenome
- results/analysis/summary_tables
- results/plots/downstream
