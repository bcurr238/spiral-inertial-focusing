// Headless copy of the page's particle model (index.html step()), run on Martel & Toner (2012) devices.
// Geometry assumptions (not given in the paper text): 60 mm long out-spiral, R grows linearly 3 -> 7 mm.
const cfg = Object.assign({ H: 0.05, Ws: [0.05, 0.075, 0.1, 0.2, 0.3, 0.4], Us: [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75],
  aUm: 10, L: 60, Rs: 3, Re_: 7, n: 500, D: 1e-4, stations: 30, dir: 'out' }, JSON.parse(arguments[0] || '{}'));
const RHO = 1000, MU = 1e-3, CLA = 0.5, RE_REF = 50, ZEQ = 0.6, YW = 0.6, CY = 3, CZ = 1, LZN = 1 / (0.385 * ZEQ ** 3);
const coreOf = (k) => (k > 1e-6 ? (-2 + Math.sqrt(4 + k * k)) / k : 0);
function run(W, Ums) {
  const H = cfg.H, DH = 2 * W * H / (W + H), HW = H / W, U = Ums, Umm = U * 1000;
  const Re = RHO * U * DH * 1e-3 / MU, CL = CLA * Math.pow(Re / RE_REF, -0.5);
  const a = cfg.aUm * 1e-6, ul = 1000 * CL * RHO * U * U * a ** 3 / (3 * Math.PI * MU * (DH * 1e-3) ** 2);
  const Rof = (s) => cfg.dir === 'out' ? cfg.Rs + (cfg.Re_ - cfg.Rs) * s / cfg.L : cfg.Re_ - (cfg.Re_ - cfg.Rs) * s / cfg.L;
  const ymax = 1 - cfg.aUm * 1e-3 / W, zmax = 1 - cfg.aUm * 1e-3 / H;
  const rate = ul * CZ * LZN * 2 * ZEQ * ZEQ * 2 / H + 1e-9;
  const h = Math.min(1.5e-4, 0.4 / rate, 0.02 * W / (Umm + 1));
  const sig = Math.sqrt(2 * cfg.D * h) * 3.4641;
  const NS = cfg.stations, snaps = Array.from({ length: NS }, () => []);
  for (let i = 0; i < cfg.n; i++) {
    let Y = -0.97 + 1.94 * Math.random(), Z = (2 * Math.random() - 1) * 0.88, s = 0, k = 0, guard = 0;
    Y = Math.max(-ymax, Math.min(ymax, Y)); Z = Math.max(-zmax, Math.min(zmax, Z));
    while (s < cfg.L && guard++ < 5e6) {
      const R = Rof(s), De = Re * Math.sqrt(DH / (2 * R)), kap = Math.min(1.2, 0.035 * De), yc = coreOf(kap);
      const UD = 0.18 * De ** 1.63, UDN = UD / ((1 - yc * yc) ** 2 * Math.exp(kap * yc)), G = Math.min(0.3, 0.01 * De);
      const Y2 = Y * Y, Z2 = Z * Z, oy = 1 - Y2, oz = 1 - Z2, e = UDN * Math.exp(kap * Y);
      let vy = e * oy * oy * oz * (1 - 5 * Z2), vz = -e * HW * oy * (kap * oy - 4 * Y) * Z * oz * oz;
      const t = Math.abs(Y); if (t > YW) { const p = (t - YW) / (1 - YW); vy -= Math.sign(Y) * ul * CY * p * p; }
      vz -= ul * CZ * Z * (Z2 - ZEQ * ZEQ) * LZN;
      Y += vy * h * 2 / W + (Math.random() - 0.5) * sig * 2 / W; Z += vz * h * 2 / H + (Math.random() - 0.5) * sig * 2 / H;
      if (Y > ymax) Y = ymax; else if (Y < -ymax) Y = -ymax; if (Z > zmax) Z = zmax; else if (Z < -zmax) Z = -zmax;
      const u = Math.max(0.04, 1.6875 * (1 - Y ** 8) * (1 + G * Y) * (1 - Z * Z)) * Umm;
      s += u * h * R / (R + Y * W / 2);
      while (k < NS && s >= (k + 1) * cfg.L / NS) { snaps[k].push([(Y + 1) / 2 * W * 1000, Z]); k++; }
    }
  }
  return { Re, ReC: 1.5 * Re, DeMid: Re * Math.sqrt(DH / (cfg.Rs + cfg.Re_)), ul, UDmid: 0.18 * (Re * Math.sqrt(DH / (cfg.Rs + cfg.Re_))) ** 1.63, snaps };
}
// Streak metrics in the paper's style: lateral intensity profile, 5% threshold = width, 50% region centre = position
function streak(pts, Wum) {
  const bw = 2, nb = Math.ceil(Wum / bw), hist = new Float64Array(nb);
  for (const p of pts) hist[Math.min(nb - 1, Math.floor(p[0] / bw))]++;
  const sm = hist.map((v, i) => (0.25 * (hist[i - 1] || 0) + 0.5 * v + 0.25 * (hist[i + 1] || 0)));
  let mx = 0, im = 0; sm.forEach((v, i) => { if (v > mx) { mx = v; im = i; } });
  let lo = nb, hi = -1; sm.forEach((v, i) => { if (v >= 0.05 * mx) { lo = Math.min(lo, i); hi = Math.max(hi, i); } });
  let l5 = im, h5 = im; while (l5 > 0 && sm[l5 - 1] >= 0.5 * mx) l5--; while (h5 < nb - 1 && sm[h5 + 1] >= 0.5 * mx) h5++;
  let peaks = 0, above = false; sm.forEach((v) => { const a = v >= 0.25 * mx; if (a && !above) peaks++; above = a; });
  const zs = pts.map((p) => Math.abs(p[1])); zs.sort((x, y) => x - y);
  return { width: (hi - lo + 1) * bw, pos: ((l5 + h5 + 1) / 2) * bw, peaks, absZmed: zs[zs.length >> 1] };
}
const out = [];
for (const W of cfg.Ws) for (const U of cfg.Us) {
  const r = run(W, U), Wum = W * 1000, a = cfg.aUm;
  const st = r.snaps.map((pts) => streak(pts, Wum)), end = st[st.length - 1];
  let Lf = null; for (let k = 0; k < st.length; k++) if (st.slice(k).every((q) => q.width < 3.8 * a)) { Lf = (k + 1) * cfg.L / st.length; break; }
  out.push({ W: Wum, U, Re: +r.Re.toFixed(0), ReC: +r.ReC.toFixed(0), De: +r.DeMid.toFixed(1), liftDean: +(r.ul / r.UDmid).toFixed(2),
    width_a: +(end.width / a).toFixed(1), pos_a: +(end.pos / a).toFixed(1), peaks: end.peaks, focused: end.width < 3.8 * a, Lfocus_mm: Lf, absZ: +end.absZmed.toFixed(2) });
}
print(JSON.stringify(out));
