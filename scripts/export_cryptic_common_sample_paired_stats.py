import csv
import sys
from collections import defaultdict
from itertools import combinations
from math import comb

RPT = "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication"
TBP = "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据\\bdq_calls.csv"
UP = {"p.Ile67Leu", "p.Phe19Ser", "p.Ser68Asn"}
DN = {"p.Glu55Asp"}


def load_variant_calls():
    by = defaultdict(list)
    with open(
        RPT + "\\final_external_bdq_phenotype_validation_variant_calls.csv",
        encoding="utf-8-sig",
    ) as f:
        for r in csv.DictReader(f):
            p = (r.get("amr_hunter_prediction") or "").strip().upper()
            if p not in ("R", "S"):
                p = "UNKNOWN"
            g = (r.get("gene") or "").strip()
            m = (r.get("mutation") or "").strip()
            vc = (r.get("variant_class") or "").strip()
            ev = (r.get("evidence_code") or "").strip()
            wid = (r.get("matched_who_row_id") or "").strip()
            if g == "Rv0678" and m in UP:
                p = "R"
            if g == "Rv0678" and m in DN and (p == "R"):
                p = "UNKNOWN"
            if g == "mtrB" and ev == "MTRAB_DEREGULATION":
                p = "S"
            if (
                g == "pepQ"
                and ev == "NEGATIVE_LOF_REVIEWED"
                and (vc == "protein_substitution")
            ):
                p = "UNKNOWN"
            by[r["unique_id"]].append((g, m, vc, ev, p, wid))
    return by


def load_pheno():
    pheno = {}
    with open(
        RPT + "\\final_external_bdq_phenotype_validation_isolates.csv",
        encoding="utf-8-sig",
    ) as f:
        for r in csv.DictReader(f):
            pheno[r["unique_id"]] = (r["bdq_binary_phenotype"], r.get("bdq_mic") or "")
    return pheno


def load_tbp():
    calls = {}
    with open(TBP, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            calls[r["id"]] = (r["tbprofiler_bdq"] or "").strip()
    return calls


def load_catalogue():
    grade = {}
    with open(RPT + "\\who_bdq_per_row_benchmark.csv", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            grade[r["who_row_id"]] = r["who_truth"]
    iso_cat = defaultdict(set)
    with open(
        RPT + "\\final_external_bdq_phenotype_validation_variant_calls.csv",
        encoding="utf-8-sig",
    ) as f:
        for r in csv.DictReader(f):
            wid = (r.get("matched_who_row_id") or "").strip()
            if wid and wid in grade:
                iso_cat[r["unique_id"]].add(grade[wid])
    return {u: "R" if "R" in s else "S" for u, s in iso_cat.items()}


def amr_isolate_calls(by):
    iso = {}
    for u, rows in by.items():
        ps = [x[4] for x in rows]
        iso[u] = "R" if "R" in ps else "S" if "S" in ps else "UNKNOWN"
    return iso


def exact_mcnemar(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return 2 * sum((comb(n, i) for i in range(k + 1))) * 0.5**n


def holm(ps):
    m = len(ps)
    idx = sorted(range(m), key=lambda i: ps[i])
    adj = [1.0] * m
    prev = 0.0
    for rank, i in enumerate(idx):
        raw = ps[i]
        adj[i] = min(1.0, max(prev, raw * (m - rank)))
        prev = max(prev, adj[i])
    return adj


def cochran_q(success_by_method, n_subjects):
    k = len(success_by_method)
    C = [sum(v) for v in success_by_method]
    R = [sum((v[i] for v in success_by_method)) for i in range(n_subjects)]
    R2 = sum((r * r for r in R))
    num = (k - 1) * (k * sum((c * c for c in C)) - sum(C) ** 2)
    den = k * sum(C) - R2
    return num / den if den else 0.0


by = load_variant_calls()
pheno = load_pheno()
tbp = load_tbp()
catalog = load_catalogue()
amr = amr_isolate_calls(by)
ph_ids = set(pheno)
tbp_repaired = {}
for uid, call in tbp.items():
    if uid in ph_ids:
        tbp_repaired[uid] = call
    else:
        candidates = [p for p in ph_ids if p.endswith(uid)]
        assert len(candidates) == 1, f"cannot repair TBP id {uid!r}: {candidates}"
        tbp_repaired[candidates[0]] = call
tbp = tbp_repaired
print(f"TBP ids repaired: {len(tbp_repaired)}")
tp = [u for u in pheno if pheno[u][0] == "R" and amr[u] == "R"]
fp = [u for u in pheno if pheno[u][0] == "S" and amr[u] == "R"]
fn = [u for u in pheno if pheno[u][0] == "R" and amr[u] != "R"]
un = [u for u in pheno if amr[u] == "UNKNOWN"]
assert len(tp) == 24 and len(fp) == 81 and (len(fn) == 47), (len(tp), len(fp), len(fn))
print(f"AMR full cohort: TP={len(tp)} FP={len(fp)} FN={len(fn)} UNKNOWN={len(un)}")
common = sorted(set(pheno) & set(tbp))
print(f"common isolates: {len(common)} (TB-Profiler output available)")
flags = {"AMR": {}, "TBP": {}, "CAT": {}}
for u in common:
    flags["AMR"][u] = 1 if amr[u] == "R" else 0
    flags["TBP"][u] = 1 if tbp[u] == "R" else 0
    flags["CAT"][u] = 1 if catalog.get(u) == "R" else 0
methods = ["AMR", "TBP", "CAT"]
r_sub = [u for u in common if pheno[u][0] == "R"]
s_sub = [u for u in common if pheno[u][0] == "S"]
print(f"common R-phenotype: {len(r_sub)}; common S-phenotype: {len(s_sub)}")
conf = {m: [0, 0, 0, 0] for m in methods}
for m in methods:
    for u in common:
        ph = pheno[u][0]
        cc = flags[m][u]
        if ph == "R":
            conf[m][0 if cc else 3] += 1
        else:
            conf[m][1 if not cc else 2] += 1
rows = []
for m in methods:
    tp_, tn_, fp_, fn_ = conf[m]
    sens = 100 * tp_ / (tp_ + fn_)
    spec = 100 * tn_ / (tn_ + fp_)
    rows.append([f"common_{len(common)}", m, tp_, tn_, fp_, fn_, sens, spec])
    print(
        f"{m:>4}: TP={tp_:>3} TN={tn_:>4} FP={fp_:>3} FN={fn_:>2} sens={sens:.1f}% spec={spec:.3f}%"
    )
q_sens = cochran_q([[flags[m][u] for u in r_sub] for m in methods], len(r_sub))
q_spec = cochran_q([[1 - flags[m][u] for u in s_sub] for m in methods], len(s_sub))
print(f"Cochran Q (sensitivity, common {len(common)}): Q={q_sens:.2f}, df=2")
print(f"Cochran Q (specificity, common {len(common)}): Q={q_spec:.2f}, df=2")
labels = {}
for m1, m2 in combinations(methods, 2):
    b_s = c_s = b_p = c_p = 0
    for u in r_sub:
        if flags[m1][u] and (not flags[m2][u]):
            b_s += 1
        elif not flags[m1][u] and flags[m2][u]:
            c_s += 1
    for u in s_sub:
        if flags[m1][u] and (not flags[m2][u]):
            b_p += 1
        elif not flags[m1][u] and flags[m2][u]:
            c_p += 1
    p_sens = exact_mcnemar(b_s, c_s)
    p_spec = exact_mcnemar(b_p, c_p)
    labels[m1, m2] = (b_s, c_s, p_sens, b_p, c_p, p_spec)
    print(
        f"{m1} vs {m2}: sens discord b={b_s} c={c_s} p={p_sens:.4g} | spec discord b={b_p} c={c_p} p={p_spec:.4g}"
    )
raw_sens = [labels[p][2] for p in labels]
raw_spec = [labels[p][5] for p in labels]
adj_sens = dict(zip(labels.keys(), holm(raw_sens)))
adj_spec = dict(zip(labels.keys(), holm(raw_spec)))
print(
    "\nHolm-adjusted pairwise P (sensitivity):",
    {f"{a}-{b}": round(v, 4) for (a, b), v in adj_sens.items()},
)
print(
    "Holm-adjusted pairwise P (specificity):",
    {f"{a}-{b}": round(v, 4) for (a, b), v in adj_spec.items()},
)
with open(
    RPT + "\\cryptic_common_sample_paired_stats.csv",
    "w",
    newline="",
    encoding="utf-8-sig",
) as f:
    w = csv.writer(f)
    w.writerow(["analysis", "method_a", "method_b", "statistic", "value"])
    w.writerow(
        [
            f"common_{len(common)}",
            "all",
            "all",
            "cochran_q_sensitivity",
            f"{q_sens:.4f}",
        ]
    )
    w.writerow(
        [
            f"common_{len(common)}",
            "all",
            "all",
            "cochran_q_specificity",
            f"{q_spec:.4f}",
        ]
    )
    for (a, b), v in adj_sens.items():
        w.writerow([f"common_{len(common)}", a, b, "mcnemar_sens_holm", f"{v:.6f}"])
    for (a, b), v in adj_spec.items():
        w.writerow([f"common_{len(common)}", a, b, "mcnemar_spec_holm", f"{v:.6f}"])
    for m in methods:
        w.writerow([f"common_{len(common)}", m, "", "tp", conf[m][0]])
        w.writerow([f"common_{len(common)}", m, "", "tn", conf[m][1]])
        w.writerow([f"common_{len(common)}", m, "", "fp", conf[m][2]])
        w.writerow([f"common_{len(common)}", m, "", "fn", conf[m][3]])
with open(
    RPT + "\\cryptic_common_sample_per_isolate_calls.csv",
    "w",
    newline="",
    encoding="utf-8-sig",
) as f:
    w = csv.writer(f)
    w.writerow(["unique_id", "phenotype", "amr_call", "tbp_call", "catalogue_call"])
    for u in common:
        w.writerow([u, pheno[u][0], amr[u], tbp[u], catalog.get(u, "")])
with open(
    RPT + "\\cryptic_final_per_isolate_calls_20260914.csv",
    "w",
    newline="",
    encoding="utf-8-sig",
) as f:
    w = csv.writer(f)
    w.writerow(
        [
            "unique_id",
            "bdq_binary_phenotype",
            "bdq_mic",
            "amr_final_call",
            "tbp_call",
            "catalogue_call",
        ]
    )
    for u in sorted(pheno):
        w.writerow(
            [u, pheno[u][0], pheno[u][1], amr[u], tbp.get(u, ""), catalog.get(u, "")]
        )
print(
    "saved: cryptic_common_sample_paired_stats.csv / cryptic_common_sample_per_isolate_calls.csv"
)
print("saved: cryptic_final_per_isolate_calls_20260914.csv (full 8,535 cohort)")
