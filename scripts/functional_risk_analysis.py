#!/usr/bin/env python3
"""
functional_risk_analysis.py
Phân tích mức độ nguy hiểm của lỗi FP/FN từ 4 variant caller.

Hai lớp đánh giá:
  Lớp 1 — SnpEff Impact Distribution: phân loại lỗi theo HIGH/MODERATE/LOW/MODIFIER
  Lớp 2 — ACMG 5-tier Classification: phân loại lâm sàng P/LP/VUS/LB/B

Sử dụng:
    python3 scripts/functional_risk_analysis.py

Input:  results/annotation/{caller}/{fn|fp}.annotated.vcf.gz
Output:
    results/annotation/impact_distribution.csv    — phân bố impact per caller
    results/annotation/acmg_summary.csv           — thống kê ACMG per caller
    results/annotation/{caller}/{fn|fp}.acmg.tsv  — chi tiết từng biến thể
"""

import gzip
import csv
import os
from collections import defaultdict

# ─── Cấu hình ────────────────────────────────────────────────────────────────

ANNOT_DIR = "results/annotation"
CALLERS = ["gatk", "deepvariant", "strelka2", "freebayes"]
CATEGORIES = ["fn", "fp"]

# Ngưỡng đánh giá chức năng
SIFT_CUTOFF = 0.05        # ≤ 0.05 = Damaging
PP2_CUTOFF = 0.85         # ≥ 0.85 = Probably damaging
PROVEAN_CUTOFF = -2.5     # ≤ -2.5 = Deleterious
CADD_CUTOFF = 20          # ≥ 20 = Top 1% deleterious
METALR_CUTOFF = 0.5       # ≥ 0.5 = Damaging
MT_DAMAGING = {"D", "A"}  # D: disease_causing, A: disease_causing_automatic


# ─── Hàm phụ trợ ─────────────────────────────────────────────────────────────

def parse_info(info_str):
    """Parse chuỗi INFO thành dict."""
    d = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def parse_ann(ann_str):
    """Parse trường ANN của SnpEff, trả về dict annotation đầu tiên."""
    fields = ann_str.split("|")
    return {
        "allele":     fields[0] if len(fields) > 0 else "",
        "effect":     fields[1] if len(fields) > 1 else "",
        "impact":     fields[2] if len(fields) > 2 else "",
        "gene":       fields[3] if len(fields) > 3 else "",
        "gene_id":    fields[4] if len(fields) > 4 else "",
        "feature":    fields[5] if len(fields) > 5 else "",
        "feature_id": fields[6] if len(fields) > 6 else "",
        "biotype":    fields[7] if len(fields) > 7 else "",
        "hgvs_c":     fields[9] if len(fields) > 9 else "",
        "hgvs_p":     fields[10] if len(fields) > 10 else "",
    }


def safe_float(val):
    """Chuyển chuỗi sang float, trả về None nếu rỗng/NA."""
    try:
        val = str(val).replace(";", ",").split(",")[0].strip()
        if val in (".", "", "NA", "None"):
            return None
        return float(val)
    except (ValueError, TypeError):
        return None


def evaluate_functional(info):
    """Đánh giá chức năng bằng 6 công cụ, trả về dict tool -> (score, pred)."""
    results = {}

    # SIFT
    s = safe_float(info.get("dbNSFP_SIFT_score", "."))
    results["SIFT"] = (s, "D" if s is not None and s <= SIFT_CUTOFF else
                          (info.get("dbNSFP_SIFT_pred", ".") if s is None else "T"))

    # PolyPhen-2
    s = safe_float(info.get("dbNSFP_Polyphen2_HDIV_score", "."))
    results["PolyPhen2"] = (s, "D" if s is not None and s >= PP2_CUTOFF else
                              (info.get("dbNSFP_Polyphen2_HDIV_pred", ".") if s is None else "T"))

    # PROVEAN
    s = safe_float(info.get("dbNSFP_PROVEAN_score", "."))
    results["PROVEAN"] = (s, "D" if s is not None and s <= PROVEAN_CUTOFF else
                            (info.get("dbNSFP_PROVEAN_pred", ".") if s is None else "N"))

    # MutationTaster
    mt_score = safe_float(info.get("dbNSFP_MutationTaster_score", "."))
    mt_pred_raw = info.get("dbNSFP_MutationTaster_pred", ".")
    mt_pred = mt_pred_raw.split(",")[0].strip() if mt_pred_raw not in (".", "") else None
    results["MutationTaster"] = (mt_score, mt_pred)

    # CADD
    s = safe_float(info.get("dbNSFP_CADD_phred", "."))
    results["CADD"] = (s, "D" if s is not None and s >= CADD_CUTOFF else
                         (None if s is None else "T"))

    # MetaLR
    s = safe_float(info.get("dbNSFP_MetaLR_score", "."))
    results["MetaLR"] = (s, "D" if s is not None and s >= METALR_CUTOFF else
                           (info.get("dbNSFP_MetaLR_pred", ".") if s is None else "T"))

    return results


def count_damaging(func_results):
    """Đếm số công cụ dự đoán biến thể có hại."""
    count = 0
    for tool, (score, pred) in func_results.items():
        if tool == "MutationTaster":
            if pred in MT_DAMAGING:
                count += 1
        elif pred == "D":
            count += 1
    return count


def classify_acmg(info, ann, func_results):
    """
    Phân loại ACMG 5-tier (đơn giản hoá).

    Evidence:
      PVS1 — null variant (nonsense, frameshift, splice) / HIGH impact
      PM2  — biến thể hiếm (MAF < 1% hoặc absent)
      PP3  — ≥3 công cụ in-silico dự đoán có hại
      PP5  — ClinVar Pathogenic
      BP4  — tất cả công cụ đều dự đoán lành tính
    """
    ev_p, ev_b = [], []

    impact = ann.get("impact", "")
    effect = ann.get("effect", "")
    clnsig = info.get("CLNSIG", "")
    maf = safe_float(info.get("dbNSFP_1000Gp3_AF", "."))

    n_damaging = count_damaging(func_results)
    n_with_data = sum(1 for _, (s, p) in func_results.items() if s is not None or p is not None)

    # PVS1: null variant / HIGH impact
    null_effects = {"stop_gained", "frameshift_variant", "splice_acceptor_variant",
                    "splice_donor_variant", "start_lost", "stop_lost"}
    if effect in null_effects or impact == "HIGH":
        ev_p.append("PVS1")

    # PM2: hiếm hoặc absent
    if maf is None or maf < 0.01:
        ev_p.append("PM2")

    # PP3: ≥3 công cụ dự đoán có hại
    if n_damaging >= 3:
        ev_p.append("PP3")

    # BP4: tất cả lành tính
    if n_with_data >= 3 and n_damaging == 0:
        ev_b.append("BP4")

    # PP5: ClinVar
    if "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig:
        ev_p.append("PP5")

    # --- Quy tắc phân loại ---
    has = lambda tag, lst: tag in lst
    pvs1 = has("PVS1", ev_p)
    pm2  = has("PM2", ev_p)
    pp3  = has("PP3", ev_p)
    pp5  = has("PP5", ev_p)
    bp4  = has("BP4", ev_b)

    if pvs1 and (pm2 or pp3 or pp5):
        acmg = "P"
    elif pvs1:
        acmg = "LP"
    elif pm2 and pp3:
        acmg = "LP"
    elif bp4 and not pvs1 and not pp3 and maf is not None and maf >= 0.01:
        acmg = "B"
    elif bp4 and not pvs1 and not pp3:
        acmg = "LB"
    else:
        acmg = "VUS"

    return acmg, ev_p, ev_b


# ─── Xử lý chính ─────────────────────────────────────────────────────────────

def process_vcf(vcf_path):
    """Đọc VCF annotated, trả về danh sách record đã phân loại."""
    records = []
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            info = parse_info(cols[7])

            # Parse ANN
            ann_raw = info.get("ANN", "")
            ann = parse_ann(ann_raw.split(",")[0]) if ann_raw else {}

            # Đánh giá chức năng
            func = evaluate_functional(info)
            n_damaging = count_damaging(func)

            # Phân loại ACMG
            acmg_class, ev_p, ev_b = classify_acmg(info, ann, func)

            records.append({
                "CHROM": chrom,
                "POS": pos,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "GENE": ann.get("gene", ""),
                "EFFECT": ann.get("effect", ""),
                "IMPACT": ann.get("impact", ""),
                "HGVS_C": ann.get("hgvs_c", ""),
                "HGVS_P": ann.get("hgvs_p", ""),
                "SIFT": func["SIFT"][0],
                "SIFT_PRED": func["SIFT"][1],
                "PP2": func["PolyPhen2"][0],
                "PP2_PRED": func["PolyPhen2"][1],
                "PROVEAN": func["PROVEAN"][0],
                "PROVEAN_PRED": func["PROVEAN"][1],
                "MT": func["MutationTaster"][0],
                "MT_PRED": func["MutationTaster"][1],
                "CADD": func["CADD"][0],
                "CADD_PRED": func["CADD"][1],
                "METALR": func["MetaLR"][0],
                "METALR_PRED": func["MetaLR"][1],
                "MAF_1KG": safe_float(info.get("dbNSFP_1000Gp3_AF", ".")),
                "CLINVAR": info.get("CLNSIG", "."),
                "N_DAMAGING": n_damaging,
                "ACMG_EV_P": ";".join(ev_p) if ev_p else ".",
                "ACMG_EV_B": ";".join(ev_b) if ev_b else ".",
                "ACMG": acmg_class,
            })

    return records


def write_tsv(records, output_path):
    """Ghi danh sách record ra file TSV."""
    if not records:
        print(f"  (trống) -> bỏ qua {output_path}")
        return

    fieldnames = list(records[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            row = {k: (v if v is not None else ".") for k, v in rec.items()}
            writer.writerow(row)

    print(f"  -> {output_path}  ({len(records)} variants)")


def main():
    impact_rows = []
    acmg_rows = []

    for caller in CALLERS:
        for cat in CATEGORIES:
            vcf_path = os.path.join(ANNOT_DIR, caller, f"{cat}.annotated.vcf.gz")

            if not os.path.exists(vcf_path):
                print(f"  Bỏ qua {caller}/{cat} — chưa có file annotated")
                continue

            print(f"\n{'='*55}")
            print(f"  {caller} / {cat.upper()}")
            print(f"  Input: {vcf_path}")
            print(f"{'='*55}")

            records = process_vcf(vcf_path)

            # --- Ghi chi tiết ---
            out_tsv = os.path.join(ANNOT_DIR, caller, f"{cat}.acmg.tsv")
            write_tsv(records, out_tsv)

            # --- Lớp 1: SnpEff Impact Distribution ---
            impact_counts = defaultdict(int)
            for rec in records:
                impact_counts[rec["IMPACT"]] += 1

            impact_rows.append({
                "Caller": caller,
                "Category": cat.upper(),
                "Total": len(records),
                "HIGH": impact_counts.get("HIGH", 0),
                "MODERATE": impact_counts.get("MODERATE", 0),
                "LOW": impact_counts.get("LOW", 0),
                "MODIFIER": impact_counts.get("MODIFIER", 0),
            })

            # --- Lớp 2: ACMG Classification ---
            acmg_counts = defaultdict(int)
            for rec in records:
                acmg_counts[rec["ACMG"]] += 1

            acmg_rows.append({
                "Caller": caller,
                "Category": cat.upper(),
                "Total": len(records),
                "P": acmg_counts.get("P", 0),
                "LP": acmg_counts.get("LP", 0),
                "VUS": acmg_counts.get("VUS", 0),
                "LB": acmg_counts.get("LB", 0),
                "B": acmg_counts.get("B", 0),
            })

    # ─── Ghi Impact Distribution CSV ─────────────────────────────────────
    impact_path = os.path.join(ANNOT_DIR, "impact_distribution.csv")
    with open(impact_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Caller", "Category", "Total", "HIGH", "MODERATE", "LOW", "MODIFIER"])
        writer.writeheader()
        writer.writerows(impact_rows)

    # ─── Ghi ACMG Summary CSV ────────────────────────────────────────────
    acmg_path = os.path.join(ANNOT_DIR, "acmg_summary.csv")
    with open(acmg_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Caller", "Category", "Total", "P", "LP", "VUS", "LB", "B"])
        writer.writeheader()
        writer.writerows(acmg_rows)

    # ─── In kết quả ──────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print("  LỚP 1: SnpEff Impact Distribution")
    print(f"{'='*65}")
    print(f"{'Caller':<15} {'Cat':<4} {'Total':>6} {'HIGH':>6} {'MOD':>6} {'LOW':>6} {'MODIF':>6}")
    print("-" * 65)
    for r in impact_rows:
        print(f"{r['Caller']:<15} {r['Category']:<4} {r['Total']:>6} "
              f"{r['HIGH']:>6} {r['MODERATE']:>6} {r['LOW']:>6} {r['MODIFIER']:>6}")

    print(f"\n{'='*65}")
    print("  LỚP 2: ACMG Classification")
    print(f"{'='*65}")
    print(f"{'Caller':<15} {'Cat':<4} {'Total':>6} {'P':>4} {'LP':>4} {'VUS':>5} {'LB':>4} {'B':>4}")
    print("-" * 65)
    for r in acmg_rows:
        print(f"{r['Caller']:<15} {r['Category']:<4} {r['Total']:>6} "
              f"{r['P']:>4} {r['LP']:>4} {r['VUS']:>5} {r['LB']:>4} {r['B']:>4}")

    print(f"\n  Output files:")
    print(f"    {impact_path}")
    print(f"    {acmg_path}")
    print(f"    results/annotation/{{caller}}/{{fn|fp}}.acmg.tsv")


if __name__ == "__main__":
    main()
