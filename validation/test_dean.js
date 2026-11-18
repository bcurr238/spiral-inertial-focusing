load('flow_core.js');
for (const [W, H, NY, NZ] of [[0.5, 0.13, 96, 24], [0.5, 0.13, 128, 32], [2.0, 0.13, 256, 24], [0.4, 0.05, 128, 16], [0.3, 0.2, 48, 32], [0.8, 0.06, 160, 16]]) {
  const t0 = Date.now(), F = solveDean(W, H, NY, NZ), o = {};
  sampleDean(F, 0, 0, o);
  const wide = 0.0238 * (H / 2) ** 2 * 2.25;   // analytic centre value for an infinitely wide channel (u_max = 1.5 U)
  let zc = 0; for (let z = 0.05; z < 1; z += 0.005) { sampleDean(F, 0, z, o); if (o.vy < 0) { zc = z; break; } }
  sampleDean(F, 0, 0, o); const c = o.vy; sampleDean(F, 0.5, 0, o); const h = o.vy; sampleDean(F, 0.8, 0, o); const e = o.vy; sampleDean(F,0.5,0,o); const sh = o.u1;
  print(`${W*1000}x${H*1000} grid ${NY}x${NZ}: ${Date.now()-t0} ms, ${F.iters} its | vy(0,0)=${c.toExponential(3)} wide-limit=${wide.toExponential(3)} ratio=${(c/wide).toFixed(3)} | zero-cross Z=${zc.toFixed(3)} | vy at Y=0.5: ${(h/c).toFixed(2)} Y=0.8: ${(e/c).toFixed(2)} of centre | core Y=${F.coreY.toFixed(2)} Z=${F.coreZ.toFixed(2)} | uMax/U=${F.uMax.toFixed(3)} vAvg=${F.vAvg.toExponential(3)} vAvg/vOut=${(F.vAvg/F.vOut).toFixed(2)} | u1(0.5,0)=${sh.toExponential(2)}`);
}
