// ── Flow core: Dean flow solver and lift field (shared verbatim by the page and the validation harness) ──
// Leading-order secondary flow in a curved rectangular duct, valid at low Dean number:
//   primary flow      ∇²u0 = −G/μ                          (straight-duct Poiseuille flow)
//   secondary flow    μ∇⁴ψ = −(ρ/R) ∂(u0²)/∂z ,  ψ = ∂ψ/∂n = 0 on the walls,  v_y = ∂ψ/∂z, v_z = −∂ψ/∂y
//   primary shift     μ∇²u1 = ρ (v·∇)u0                    (outward shift of the velocity peak)
// Everything is solved once per cross-section for unit ρU²/(μR); the fields then scale linearly with it.
// Lengths in mm. Returned fields live on an (NY+1)×(NZ+1) node grid, index i*(NZ+1)+j, Y and Z from −1 to 1.
function solveDean(W, H, NY, NZ) {
  const dy = W / NY, dz = H / NZ, n1 = NZ + 1, n = (NY + 1) * n1, id = (i, j) => i * n1 + j;
  const iy2 = 1 / (dy * dy), iz2 = 1 / (dz * dz);
  const cg = (apply, b, tol, maxIt) => {            // conjugate gradients on interior nodes
    const x = new Float64Array(n), r = Float64Array.from(b), p = Float64Array.from(b), Ap = new Float64Array(n);
    let rr = 0, bb = 0; for (let k = 0; k < n; k++) { rr += r[k] * r[k]; } bb = rr || 1;
    let it = 0;
    for (; it < maxIt && rr > tol * tol * bb; it++) {
      apply(p, Ap); let pAp = 0; for (let k = 0; k < n; k++) pAp += p[k] * Ap[k];
      const al = rr / pAp; let rn = 0;
      for (let k = 0; k < n; k++) { x[k] += al * p[k]; r[k] -= al * Ap[k]; rn += r[k] * r[k]; }
      const be = rn / rr; rr = rn; for (let k = 0; k < n; k++) p[k] = r[k] + be * p[k];
    }
    x.iters = it; return x;
  };
  const negLap = (x, out) => {                       // −∇² with x = 0 on the walls
    out.fill(0);
    for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) { const k = id(i, j);
      out[k] = (2 * x[k] - x[k - n1] - x[k + n1]) * iy2 + (2 * x[k] - x[k - 1] - x[k + 1]) * iz2; }
  };
  const biharm = (x, out) => {                       // ∇⁴ with clamped walls (mirror ghost nodes)
    out.fill(0);
    for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) { const k = id(i, j), c = x[k];
      const ym2 = i === 1 ? c : x[k - 2 * n1], yp2 = i === NY - 1 ? c : x[k + 2 * n1];
      const zm2 = j === 1 ? c : x[k - 2], zp2 = j === NZ - 1 ? c : x[k + 2];
      out[k] = (ym2 - 4 * x[k - n1] + 6 * c - 4 * x[k + n1] + yp2) * iy2 * iy2
             + (zm2 - 4 * x[k - 1] + 6 * c - 4 * x[k + 1] + zp2) * iz2 * iz2
             + 2 * iy2 * iz2 * (x[k - n1 - 1] + x[k - n1 + 1] + x[k + n1 - 1] + x[k + n1 + 1] - 2 * (x[k - n1] + x[k + n1] + x[k - 1] + x[k + 1]) + 4 * c); }
  };
  const interior = (f) => { const b = new Float64Array(n); for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) b[id(i, j)] = f(i, j); return b; };
  const mean = (x) => {                              // trapezoid mean over the cross-section (walls are zero)
    let s = 0; for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) s += x[id(i, j)]; return s / (NY * NZ); };
  // primary flow, normalised to unit mean
  const u0 = cg(negLap, interior(() => 1), 1e-7, 4000); { const m = mean(u0); for (let k = 0; k < n; k++) u0[k] /= m; }
  // secondary-flow stream function per unit ρU²/(μR)
  const psi = cg(biharm, interior((i, j) => -(u0[id(i, j + 1)] ** 2 - u0[id(i, j - 1)] ** 2) / (2 * dz)), 1e-6, 40000);
  const vy = new Float64Array(n), vz = new Float64Array(n);
  for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) { const k = id(i, j);
    vy[k] = (psi[k + 1] - psi[k - 1]) / (2 * dz); vz[k] = -(psi[k + n1] - psi[k - n1]) / (2 * dy); }
  // first correction to the primary flow per unit (ρ/μ)·ρU²/(μR), zero net flux
  const u1 = cg(negLap, interior((i, j) => { const k = id(i, j);
    return -(vy[k] * (u0[k + n1] - u0[k - n1]) / (2 * dy) + vz[k] * (u0[k + 1] - u0[k - 1]) / (2 * dz)); }), 1e-7, 4000);
  { const m = mean(u1); for (let k = 0; k < n; k++) u1[k] -= m * u0[k]; }
  // local shear of the primary flow, scaled so that a parabolic profile across the height gives sz = Z
  const uc = u0[id(NY >> 1, NZ >> 1)], sy = new Float64Array(n), sz = new Float64Array(n), gY = new Float32Array(NY + 1), gZ = new Float32Array(NZ + 1);
  for (let i = 0; i <= NY; i++) for (let j = 0; j <= NZ; j++) { const k = id(i, j);
    const dudy = i === 0 ? (u0[k + n1] - u0[k]) / dy : i === NY ? (u0[k] - u0[k - n1]) / dy : (u0[k + n1] - u0[k - n1]) / (2 * dy);
    const dudz = j === 0 ? (u0[k + 1] - u0[k]) / dz : j === NZ ? (u0[k] - u0[k - 1]) / dz : (u0[k + 1] - u0[k - 1]) / (2 * dz);
    sy[k] = -(H / 2) * dudy / (2 * uc); sz[k] = -(H / 2) * dudz / (2 * uc); }
  for (let i = 0; i <= NY; i++) gY[i] = u0[id(i, NZ >> 1)] / uc;
  for (let j = 0; j <= NZ; j++) gZ[j] = u0[id(NY >> 1, j)] / uc;
  let vAvg = 0, vOut = 0, pk = 0, coreY = 0, coreZ = 0, uMax = 0;
  for (let i = 1; i < NY; i++) for (let j = 1; j < NZ; j++) { const k = id(i, j);
    vAvg += Math.hypot(vy[k], vz[k]); if (vy[k] > vOut) vOut = vy[k]; if (u0[k] > uMax) uMax = u0[k];
    if (psi[k] > pk) { pk = psi[k]; coreY = -1 + 2 * i / NY; coreZ = -1 + 2 * j / NZ; } }
  return { NY, NZ, u0: Float32Array.from(u0), u1: Float32Array.from(u1), vy: Float32Array.from(vy), vz: Float32Array.from(vz), psi: Float32Array.from(psi), sy: Float32Array.from(sy), sz: Float32Array.from(sz), gY, gZ,
    vAvg: vAvg / (NY * NZ), vOut, psiMax: pk, coreY, coreZ, uMax, iters: psi.iters };
}
// Bilinear sample of the fields at (Y, Z) in −1…1; writes o.vy, o.vz (mm² per unit ρU²/μR), o.u0, o.u1, o.sy, o.sz, o.g, o.gz
function sampleDean(F, Y, Z, o) {
  let fy = (Y + 1) * 0.5 * F.NY, fz = (Z + 1) * 0.5 * F.NZ;
  if (fy < 0) fy = 0; else if (fy > F.NY - 1e-6) fy = F.NY - 1e-6; if (fz < 0) fz = 0; else if (fz > F.NZ - 1e-6) fz = F.NZ - 1e-6;
  const i = fy | 0, j = fz | 0, a = fy - i, b = fz - j, n1 = F.NZ + 1, k = i * n1 + j;
  const w00 = (1 - a) * (1 - b), w01 = (1 - a) * b, w10 = a * (1 - b), w11 = a * b;
  o.vy = w00 * F.vy[k] + w01 * F.vy[k + 1] + w10 * F.vy[k + n1] + w11 * F.vy[k + n1 + 1];
  o.vz = w00 * F.vz[k] + w01 * F.vz[k + 1] + w10 * F.vz[k + n1] + w11 * F.vz[k + n1 + 1];
  o.u0 = w00 * F.u0[k] + w01 * F.u0[k + 1] + w10 * F.u0[k + n1] + w11 * F.u0[k + n1 + 1];
  o.u1 = w00 * F.u1[k] + w01 * F.u1[k + 1] + w10 * F.u1[k + n1] + w11 * F.u1[k + n1 + 1];
  o.sy = w00 * F.sy[k] + w01 * F.sy[k + 1] + w10 * F.sy[k + n1] + w11 * F.sy[k + n1 + 1];
  o.sz = w00 * F.sz[k] + w01 * F.sz[k + 1] + w10 * F.sz[k + n1] + w11 * F.sz[k + n1 + 1];
  o.g = F.gY[i] + a * (F.gY[i + 1] - F.gY[i]); o.gz = F.gZ[j] + b * (F.gZ[j + 1] - F.gZ[j]);
}
// Above De ≈ 17 Dean's low-curvature series stops converging and the secondary flow grows more slowly than De².
// This roll-off tends to the De^1.63 growth reported by full simulations; its exact form is an estimate.
const DE_C = 17;
function deanRollOff(De) { return Math.pow(1 + (De / DE_C) * (De / DE_C), -0.185); }
// Inertial lift, as a migration velocity. Scale (Di Carlo 2009): f_L ρ U_max² a³ / (3π μ h²), with h the channel height.
// f_L = 0.035 at Re_h = 50, falling as Re^−0.5. The value was chosen from 0.02, 0.035 and 0.07 by scoring against
// Martel & Toner (2012); it sits inside the 0.02–0.05 range usually quoted for this coefficient.
const LIFT = { FL0: 0.035, RE_L: 50, RE_EXP: -0.5, ZEQ: 0.6, ZETA0: 0.45, BETA: 0.1 };
const LZN = 1 / (0.385 * LIFT.ZEQ * LIFT.ZEQ * LIFT.ZEQ);
function liftScale(U, uMaxRatio, H, aUm) {      // mm/s; U in m/s, H in mm, a in µm
  const Um = U * uMaxRatio, ReH = 1000 * Um * H * 1e-3 / 1e-3, a = aUm * 1e-6;
  return 1000 * LIFT.FL0 * Math.pow(ReH / LIFT.RE_L, LIFT.RE_EXP) * 1000 * Um * Um * a * a * a / (3 * Math.PI * 1e-3 * H * H * 1e-6);
}
// Shape of the lift field at a sampled point d (from sampleDean), per unit liftScale. It is built to have the equilibria
// reported for straight rectangular ducts: an unstable centre, stable points near the face centres at 0.6 of the
// half-height, and saddles near the corners. Two contributions in each direction:
//   shear-gradient lift  follows the local velocity gradient toward the walls, and reverses weakly (BETA) beside a wall,
//                        which is the slow second-stage drift toward the face centres
//   wall lift            pushes away from a wall, within half a channel height of it
// Writes o.shy, o.shz (shear-gradient) and o.wy, o.wz (wall).
function liftShape(Y, Z, d, W, H, o) {
  const aY = Y < 0 ? -Y : Y, zeta = Math.max(0, 1 - (1 - aY) * W / H), zs = Y < 0 ? -zeta : zeta;
  const B = LIFT.BETA, ZEQ = LIFT.ZEQ, mz = Math.max(-B, 1 - (1 + B) * (zeta / LIFT.ZETA0) * (zeta / LIFT.ZETA0)), my = Math.max(-B, 1 - (1 + B) * (Z / ZEQ) * (Z / ZEQ));
  o.shz = LZN * ZEQ * ZEQ * d.g * d.sz * mz; o.wz = -LZN * d.g * d.g * Z * Z * Z;
  o.shy = LZN * ZEQ * ZEQ * d.gz * d.sy * my; o.wy = -LZN * d.gz * d.gz * zs * zs * zs;
}
// ── end flow core ─────────────────────────────────────────────────────────────────────────────────────
