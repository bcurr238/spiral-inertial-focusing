import json, sys
d = json.load(open(sys.argv[1])); print(sys.argv[2] if len(sys.argv) > 2 else '')
if 'W' in d[0]:
    Ws = sorted({r['W'] for r in d}); print('U m/s  ' + ''.join(f'{w:>11.0f}µm' for w in Ws) + '   (streak width | position from inner wall, in bead diameters; F = focused)')
    for U in sorted({r['U'] for r in d}):
        print(f'{U:5.2f}  ' + ''.join((lambda r: f"{r['width_a']:5.1f}|{r['pos_a']:4.1f}{' F' if r['focused'] else '  '}{('x'+str(r['peaks'])) if r.get('peaks',1)>1 else '  '} ")([x for x in d if x['W'] == w and x['U'] == U][0]) for w in Ws))
    print('median |Z| at exit, U=1 and 2.5:', [(w, [x.get('absZ') for x in d if x['W']==w and x['U'] in (1,2.5)]) for w in Ws])
    for w in Ws: print(f'{w:4.0f}µm length to focus:', [(x['U'], x['Lfocus_mm']) for x in d if x['W'] == w and x['Lfocus_mm']])
else:
    print('Q    U     Re  De   vDean   tumour      white       red         exosome    (% by outlet 1 inner … 4 outer)')
    for r in d: print(f"{r['Q']:<4} {r['U']:<5} {r['Re']:<3} {r['De']:<4} {r['vDean']:<6}  {r['ctc']:<11} {r['wbc']:<11} {r['rbc']:<11} {r['exo']:<11}")
    print('red cells in the inner quarter along the channel (16 stations):')
    for r in d: print(f"  Q={r['Q']:<4}", r['rbcInnerByLoop'])
