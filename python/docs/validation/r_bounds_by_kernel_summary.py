import csv, sys, collections, math
rows = [r for r in csv.DictReader(open(sys.argv[1]))]
err = [r for r in rows if r.get("error")]
print("fits:", len(rows), "errors:", len(err))
for r in err: print("  ERROR", r["series"], r["end"], r["flag1"], r["mode"], r["error"])
pair = collections.defaultdict(dict)
for r in rows:
    if not r.get("error"): pair[(r["series"], int(r["end"]), int(r["flag1"]))][r["mode"]] = r
names = {3: "logistic", 4: "Richards", 5: "Gompertz"}
f = lambda r, k: float(r[k])
tot = collections.Counter()
detail = "-v" in sys.argv
for key in sorted(pair, key=lambda k: (k[2], k[0], k[1])):
    p = pair[key]
    if len(p) < 2: print("unpaired", key); continue
    L, N = p["legacy"], p["new"]
    d = f(N, "aicc") - f(L, "aicc")
    k = names[key[2]]
    tot[k, "n"] += 1
    tot[k, "better"] += d < -2
    tot[k, "same"] += abs(d) <= 2
    tot[k, "worse"] += d > 2
    tot[k, "legacy_sat"] += f(L, "cum_fit_vs_data") < 0.5
    tot[k, "new_sat"] += f(N, "cum_fit_vs_data") < 0.5
    tot[k, "new_atLB"] += int(N["r_at_LB"]); tot[k, "new_atUB"] += int(N["r_at_UB"])
    tot[k, "leg_atLB"] += int(L["r_at_LB"]); tot[k, "leg_atUB"] += int(L["r_at_UB"])
    if not math.isnan(f(L, "fc_mae")):
        tot[k, "fc_n"] += 1
        tot[k, "fc_better"] += f(N, "fc_mae") < 0.98 * f(L, "fc_mae")
        tot[k, "fc_worse"] += f(N, "fc_mae") > 1.02 * f(L, "fc_mae")
    if detail or d > 2:
        print("%-8s %-20s end=%3d y0=%8.0f legacy rLB=%-8.3g | dAICc=%12.1f | r %s -> %s | K/cum %.2f -> %.2f | fcMAE %.4g -> %.4g | n %s->%s"
              % (k, key[0], key[1], f(L, "first_count"), f(L, "rLB"), d, L["r"], N["r"],
                 f(L, "cum_fit_vs_data"), f(N, "cum_fit_vs_data"), f(L, "fc_mae"), f(N, "fc_mae"), L["n"], N["n"]))
print()
for k in ("logistic", "Richards", "Gompertz"):
    if tot[k, "n"]:
        print("%-9s windows=%d  AICc lower/same/higher with new bounds: %d/%d/%d   K<0.5*cum (legacy/new): %d/%d   "
              "r on bound legacy LB/UB %d/%d, new LB/UB %d/%d   forecast MAE better/worse: %d/%d of %d"
              % (k, tot[k, "n"], tot[k, "better"], tot[k, "same"], tot[k, "worse"], tot[k, "legacy_sat"], tot[k, "new_sat"],
                 tot[k, "leg_atLB"], tot[k, "leg_atUB"], tot[k, "new_atLB"], tot[k, "new_atUB"],
                 tot[k, "fc_better"], tot[k, "fc_worse"], tot[k, "fc_n"]))
