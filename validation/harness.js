// Validation harness: the page's particle model, using the shared flow core (flow_core.js).
//   jsc harness.js -- '{"case":"martel"}'   or   '{"case":"chip","Qs":[1.75]}'   or   '{"case":"order"}'
load('flow_core.js');
const cfg = Object.assign({ case: 'martel', n: 400, dean: 'computed' }, JSON.parse(arguments[0] || '{}'));
const RHO = 1000, MU = 1e-3;
if (cfg.lift) Object.assign(LIFT, cfg.lift);
const gauss = () => (Math.random() + Math.random() + Math.random() + Math.random() - 2) * 0.866 * 2 / 1.732;
// One pass of nPart particles of diameter aUm through a spiral; returns positions at `stations` and at the exit.
function simulate(o) {
  const { W, H, U, aUm, L, Rof, n, D, init, stations, spread } = o;
  const DH = 2 * W * H / (W + H), Umm = U * 1000, Re = RHO * U * DH * 1e-3 / MU, L_ = {};
  const F = o.F, s_ = {}, NS = stations, snaps = Array.from({ length: NS }, () => []);
  for (let i = 0; i < n; i++) {
    const ai = aUm * (spread ? Math.max(0.7, Math.min(1.3, 1 + 0.1 * gauss())) : 1), a = ai * 1e-6;
    const ul = liftScale(U, o.F.uMax, H, ai);
    const ymax = 1 - ai * 1e-3 / W, zmax = 1 - ai * 1e-3 / H;
    const h = Math.min(1.5e-4, 0.4 / (ul * LZN * 3 * 2 / H + 1e-9), 0.02 * W / (Umm + 1)), sig = Math.sqrt(2 * D * h) * 3.4641;
    let [Y, Z] = init(); Y = Math.max(-ymax, Math.min(ymax, Y)); Z = Math.max(-zmax, Math.min(zmax, Z));
    let s = 0, k = 0, guard = 0;
    while (s < L && guard++ < 2e7) {
      const R = Rof(s), De = Re * Math.sqrt(DH / (2 * R)), roll = deanRollOff(De), kd = RHO * U * U / (MU * R * 1e-3);
      sampleDean(F, Y, Z, s_);
      let vy = roll * kd * 1e-3 * s_.vy, vz = roll * kd * 1e-3 * s_.vz;
      const u = Math.max(0.04, s_.u0 + roll * (RHO / MU) * kd * 1e-9 * s_.u1) * Umm;   // sampled once per sub-step, as on the page
      liftShape(Y, Z, s_, W, H, L_); vy += ul * (L_.shy + L_.wy); vz += ul * (L_.shz + L_.wz);
      Y += vy * h * 2 / W + (Math.random() - 0.5) * sig * 2 / W; Z += vz * h * 2 / H + (Math.random() - 0.5) * sig * 2 / H;
      if (Y > ymax) Y = ymax; else if (Y < -ymax) Y = -ymax; if (Z > zmax) Z = zmax; else if (Z < -zmax) Z = -zmax;
      s += u * h * R / (R + Y * W / 2);
      while (k < NS && s >= (k + 1) * L / NS) { snaps[k].push([Y, Z, ai]); k++; }
    }
  }
  return { Re, snaps, DeMid: Re * Math.sqrt(DH / (2 * Rof(L / 2))), vAvgMid: deanRollOff(Re * Math.sqrt(DH / (2 * Rof(L / 2)))) * RHO * U * U / (MU * Rof(L / 2) * 1e-3) * 1e-3 * F.vAvg };
}
// Streak as a fluorescence image would show it: each bead contributes its own projected width (a sphere's chord profile)
function streak(pts, Wum, aUm) {
  const bw = 1, nb = Math.ceil(Wum / bw), I = new Float64Array(nb);
  for (const p of pts) { const x0 = (p[0] + 1) / 2 * Wum, r = p[2] / 2;
    for (let b = Math.max(0, Math.floor((x0 - r) / bw)); b <= Math.min(nb - 1, Math.floor((x0 + r) / bw)); b++) { const dx = (b + 0.5) * bw - x0; if (Math.abs(dx) < r) I[b] += Math.sqrt(r * r - dx * dx); } }
  let mx = 0, im = 0; I.forEach((v, i) => { if (v > mx) { mx = v; im = i; } });
  let lo = nb, hi = -1; I.forEach((v, i) => { if (v >= 0.05 * mx) { lo = Math.min(lo, i); hi = Math.max(hi, i); } });
  let l5 = im, h5 = im; while (l5 > 0 && I[l5 - 1] >= 0.5 * mx) l5--; while (h5 < nb - 1 && I[h5 + 1] >= 0.5 * mx) h5++;
  let peaks = 0, up = false; I.forEach((v) => { const a = v >= 0.3 * mx; if (a && !up) peaks++; up = a; });
  const zs = pts.map((p) => Math.abs(p[1])).sort((x, y) => x - y);
  return { width: (hi - lo + 1) * bw / aUm, pos: (l5 + h5 + 1) / 2 * bw / aUm, peaks, absZ: zs[zs.length >> 1] };
}
if (cfg.case === 'order') {
  // Size ordering on the default chip with beads spread across the whole inlet (the classic 10/15/20 µm spiral test)
  const W = 0.5, H = 0.13, R0 = 3, loops = cfg.loops || 4, b = (W + 0.75) / (2 * Math.PI), th = 2 * Math.PI * loops, L = R0 * th + b * th * th / 2;
  const F = solveDean(W, H, 96, 24), rows = [];
  for (const Q of cfg.Qs || [1, 2, 3]) { const U = Q * 1.6667e-8 / (W * H * 1e-6);
    for (const a of [10, 15, 20]) {
      const r = simulate({ W, H, U, aUm: a, L, n: cfg.n, D: 1e-4, stations: 8, F, Rof: (s) => Math.sqrt(R0 * R0 + 2 * b * s), init: () => [-0.97 + 1.94 * Math.random(), (2 * Math.random() - 1) * 0.88] });
      const e = r.snaps[r.snaps.length - 1], d = e.map((p) => (p[0] + 1) / 2 * W * 1000).sort((x, y) => x - y), z = e.map((p) => Math.abs(p[1])).sort((x, y) => x - y);
      rows.push(`Q=${Q} De=${r.DeMid.toFixed(1)} a=${a}µm: distance from inner wall median ${d[d.length >> 1].toFixed(0)} µm (10–90%: ${d[d.length * 0.1 | 0].toFixed(0)}–${d[d.length * 0.9 | 0].toFixed(0)}), median |Z| ${z[z.length >> 1].toFixed(2)}`);
    } }
  print(rows.join('\n'));
} else if (cfg.case === 'martel') {
  // Martel & Toner 2012: 10 µm beads, 50 µm tall, ~60 mm long. Radii are not in the text: assumed 3 → 7 mm.
  const c = Object.assign({ H: 0.05, Ws: [0.05, 0.075, 0.1, 0.2, 0.3, 0.4], Us: [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75], L: 60, Rs: 3, Re_: 7, dir: 'out' }, cfg);
  const rows = [];
  for (const W of c.Ws) { const F = solveDean(W, c.H, Math.max(32, Math.min(160, Math.round(W / c.H * 16))), 16);
    for (const U of c.Us) {
      const r = simulate({ W, H: c.H, U, aUm: 10, L: c.L, n: c.n, D: 1e-4, stations: 30, F,
        Rof: (s) => c.dir === 'out' ? c.Rs + (c.Re_ - c.Rs) * s / c.L : c.Re_ - (c.Re_ - c.Rs) * s / c.L,
        init: () => [-0.97 + 1.94 * Math.random(), (2 * Math.random() - 1) * 0.88] });
      const st = r.snaps.map((p) => streak(p, W * 1000, 10)), e = st[st.length - 1];
      let Lf = null; for (let k = 0; k < st.length; k++) if (st.slice(k).every((q) => q.width < 3.8)) { Lf = (k + 1) * c.L / st.length; break; }
      rows.push({ W: W * 1000, U, Re: +r.Re.toFixed(0), De: +r.DeMid.toFixed(1), width_a: +e.width.toFixed(1), pos_a: +e.pos.toFixed(1), focused: e.width < 3.8, Lfocus_mm: Lf, peaks: e.peaks, absZ: +e.absZ.toFixed(2) });
    } }
  print(JSON.stringify(rows));
} else {
  // Default chip as on the page: 500 × 130 µm, 4 loops, R 3 → 8 mm, sample in the outer 15% of the width.
  const c = Object.assign({ W: 0.5, H: 0.13, R0: 3, loops: 4, Qs: [0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2, 3], YS: 0.7 }, cfg);
  const pitch = c.W + 0.75, b = pitch / (2 * Math.PI), th = 2 * Math.PI * c.loops, L = c.R0 * th + b * th * th / 2;
  const F = solveDean(c.W, c.H, 96, 24), rows = [];
  for (const Q of c.Qs) { const U = Q * 1.6667e-8 / (c.W * c.H * 1e-6), row = { Q, U: +U.toFixed(2) };
    for (const [key, a, D] of [['ctc', 20, 1e-4], ['wbc', 12, 1e-4], ['rbc', 5.6, 1e-4], ['exo', 0.1, 2e-5]].filter((x) => !cfg.species || cfg.species.includes(x[0]))) {
      const r = simulate({ W: c.W, H: c.H, U, aUm: a, L, n: c.n, D, stations: 16, F, spread: true, Rof: (s) => Math.sqrt(c.R0 * c.R0 + 2 * b * s),
        init: () => [c.YS + 0.02 + (0.96 - c.YS) * Math.random(), (2 * Math.random() - 1) * 0.88] });
      const e = r.snaps[r.snaps.length - 1], o = [0, 0, 0, 0]; e.forEach((p) => o[p[0] < -0.5 ? 0 : p[0] < 0 ? 1 : p[0] < 0.5 ? 2 : 3]++);
      row[key] = o.map((x) => Math.round(100 * x / e.length)).join('/'); row[key + '_f'] = o.map((x) => x / e.length); row.Re = +r.Re.toFixed(0); row.De = +r.DeMid.toFixed(1); row.vDean = +r.vAvgMid.toFixed(2);
      if (key === 'rbc' && !cfg.species) row.rbcInnerByLoop = r.snaps.map((p) => Math.round(100 * p.filter((q) => q[0] < -0.5).length / p.length)).join(' ');
    }
    rows.push(row); }
  print(JSON.stringify(rows));
}
