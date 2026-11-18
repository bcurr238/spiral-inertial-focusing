# Separation score for CTC isolation: collect outlet 1 (inner), or outlets 1+2.
#   recovery  = fraction of tumour cells collected
#   WBC / RBC carryover = fraction of white / red cells that end up in the collection
#   enrichment over WBC = recovery / WBC carryover (how many-fold the CTC:WBC ratio improves)
import json, sys, math
for fn in sys.argv[1:]:
    d = json.load(open(fn)); print(f'\n== {fn}')
    print('  Q     Re   De  |  collect outlet 1:  CTC rec  WBC  RBC   enrich | collect 1+2:  CTC rec  WBC   RBC   enrich')
    for r in d:
        c, w, b = r['ctc_f'], r['wbc_f'], r['rbc_f']
        def m(k):
            rec = sum(c[:k]); wb = sum(w[:k]); rb = sum(b[:k])
            en = rec / max(wb, 0.5 / 600)
            return f"{100*rec:6.0f}% {100*wb:4.0f}% {100*rb:4.0f}% {en:7.0f}x"
        print(f"  {r['Q']:<5} {r['Re']:<4} {r['De']:<5}|      {m(1)}  |     {m(2)}")
