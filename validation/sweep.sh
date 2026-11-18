#!/bin/zsh
# usage: sweep.sh OUT W H LOOPS "Q1 Q2 ..." N   — one jsc process per flow rate, run in parallel
OUT=$1; W=$2; H=$3; LOOPS=$4; QS=$5; N=$6
JSC=/System/Library/Frameworks/JavaScriptCore.framework/Versions/A/Helpers/jsc
cd "$(dirname "$0")"; mkdir -p sweep_tmp results/4-operating-sweep
for q in ${=QS}; do
  $JSC harness.js -- "{\"case\":\"chip\",\"W\":$W,\"H\":$H,\"loops\":$LOOPS,\"Qs\":[$q],\"n\":$N,\"species\":[\"ctc\",\"wbc\",\"rbc\"]}" > sweep_tmp/${OUT}_$q.json &
done
wait
python3 -I -c "
import json,glob,sys
rows=[]
for f in glob.glob('sweep_tmp/${OUT}_*.json'): rows+=json.load(open(f))
rows.sort(key=lambda r:r['Q']); json.dump(rows,open("results/4-operating-sweep/${OUT}.json","w"))"
rm -f sweep_tmp/${OUT}_*.json; rmdir sweep_tmp 2>/dev/null
