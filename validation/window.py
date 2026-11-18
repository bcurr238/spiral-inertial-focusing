# For each design: flow-rate windows where outlet 1 gives CTC recovery >= 90% with WBC and RBC carryover each <= 5%.
import json, sys
def good(r):
    c, w, b = r['ctc_f'][0], r['wbc_f'][0], r['rbc_f'][0]
    return c >= 0.90 and w <= 0.05 and b <= 0.05
for fn in sys.argv[1:]:
    d = json.load(open(fn)); wins = []; cur = None
    for r in d:
        if good(r):
            cur = cur or [r['Q'], r['Q'], []]; cur[1] = r['Q']; cur[2].append(r)
        elif cur: wins.append(cur); cur = None
    if cur: wins.append(cur)
    desc = []
    for lo, hi, rs in wins:
        best = max(rs, key=lambda r: r['ctc_f'][0] - r['wbc_f'][0] - r['rbc_f'][0])
        desc.append(f"{lo:.2f}–{hi:.2f} mL/min (width {hi-lo+0.05:.2f}, ±{100*(hi-lo+0.05)/2/((lo+hi)/2):.0f}%), best {best['Q']:.2f}: CTC {100*best['ctc_f'][0]:.0f}% WBC {100*best['wbc_f'][0]:.1f}% RBC {100*best['rbc_f'][0]:.1f}% Re {best['Re']}")
    print(f"{fn:<16} " + (' | '.join(desc) if desc else 'no window'))
