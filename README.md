# Spiral Inertial Focusing

A real-time 3D simulation of blood cells sorting by size in a spiral microchannel, driven by inertial lift and Dean flow.

**[Open the live simulation](https://bcurr238.github.io/spiral-inertial-focusing/)** (runs in the browser, no install)

This project grew out of my undergraduate and master's research in biomedical engineering: isolating circulating tumour cells (CTCs) from blood with spiral inertial focusing, and later separating exosomes from whole blood with a spiral channel combined with acoustic waves.

## What it shows

Cells enter a 500 × 130 µm, 4-loop spiral with a sheath stream and leave through four outlets. Two effects compete as they travel:

- **Inertial lift** pushes cells toward equilibrium positions near the walls. It grows with particle size as a³, so large cells focus strongly.
- **Dean flow** is a pair of counter-rotating vortices caused by the channel's curvature. Its drag sweeps small cells around the cross-section regardless of size.

Large tumour cells (20 µm) lock in near the inner wall, while red cells (5.6 µm volume-equivalent) and white cells (12 µm) ride the Dean vortices. Exosomes (100 nm) don't focus at all.

The page includes:

- four camera views (whole chip, ride the stream, a chase camera at a cross-section, the outlet fan);
- a live cross-section with the computed Dean streamlines;
- a tracked cell with force arrows for Dean drag, shear-gradient lift and wall lift, in piconewtons;
- controls for flow rate, sample concentration, inlet configuration and chip design (channel width, height, loops, inner radius);
- live readouts of Re, De, lift-to-Dean ratios and outlet percentages.

## How the physics is modelled

This is a **reduced-order educational model** that runs at 60 fps in a browser, not a CFD simulation. Each cell is tracked as a rigid sphere moving with the local flow plus its migration velocity across streamlines.

| Component | Method |
|---|---|
| Primary flow | Poiseuille flow in a rectangular duct, solved numerically on a grid |
| Dean flow | Leading-order secondary flow in a curved rectangular duct at low Dean number, solved for the actual cross-section whenever the design changes. It scales with ρU²/μR. Above De ≈ 17 an estimated roll-off slows its growth toward the De^1.63 trend from full simulations |
| Inertial lift | Magnitude f_L ρU_max²a³ / 3πμh² (Di Carlo 2009), with f_L = 0.035 at Re = 50 falling as Re^−0.5. The spatial shape is constructed to have the equilibria known for straight rectangular channels |
| Dispersion | A small constant random term standing in for cell–cell interactions |

The flow and lift code is in [`validation/flow_core.js`](validation/flow_core.js) and is pasted verbatim into `index.html`, so the page and the validation harness run identical physics.

## Validation

The model was scored against Martel & Toner (2012), who imaged 10 µm beads in 50 µm tall spirals 50–400 µm wide. It was rebuilt twice after failing this comparison.

| Their finding | Original model | Current model |
|---|---|---|
| All six widths focus 10 µm beads | Fail | Pass |
| Narrower channels focus at lower velocity | Pass | Pass |
| 400 µm channel focuses near 1.5 m/s | Fail | Pass |
| Streak moves away from the inner wall as velocity rises | Fail | Pass |
| Reversing flow direction shifts the streak toward the centre | Fail | Pass |
| Streak breaks into several near 2 m/s in wide channels | Fail | **Fail** |

On the default chip, larger beads (20 µm) focus closer to the inner wall than smaller ones (10 µm), the classic spiral result.

Raw results for each stage are in [`validation/results/`](validation/results/).

### Operating point

A sweep over flow rate, chip height and loop count found the default chip separates CTCs most robustly at **1.75 mL/min**, collecting outlet 1: about 93% tumour-cell recovery with under 1% white- or red-cell carryover, in a window of roughly 1.60–1.85 mL/min. Chips 100 µm and 160 µm tall had no clean window in the model.

## Limitations

- **The lift coefficient was calibrated on the validation data.** f_L was chosen from three candidates by scoring against Martel & Toner, so passing that comparison is not an independent test. No held-out dataset has been used yet.
- **The lift field's shape is constructed, not computed.** Its equilibria are right, but the force strength between them is interpolated.
- **No cell–cell interactions.** Sample concentration changes only how many cells are drawn, not how they focus.
- **Rigid spheres in water.** Cell deformability and blood viscosity are not modelled.
- **The operating point is a model prediction.** It has not been tested against real cell separations, and it moves if the true Dean flow strength differs.
- **High-Dean-number behaviour is approximate.** The Dean flow solution is exact only at low Dean number.

## Running it locally

The simulation is a single HTML file. Serve the folder and open it in a browser:

```bash
python3 -m http.server 8741
```

Then go to http://localhost:8741.

The validation harness runs in JavaScriptCore, which ships with macOS:

```bash
cd validation
JSC=/System/Library/Frameworks/JavaScriptCore.framework/Versions/A/Helpers/jsc
$JSC harness.js -- '{"case":"martel"}' > out.json && python3 summ.py out.json
```

Other cases are `chip` (default-chip outlets across flow rates) and `order` (size ordering). `sweep.sh`, `score.py` and `window.py` reproduce the operating-point sweep.

## References

- Martel J.M., Toner M. Inertial focusing dynamics in spiral microchannels. *Phys. Fluids* 24, 032001 (2012). https://doi.org/10.1063/1.3681228
- Valani R.N., Harding B., Stokes Y.M. Inertial particle focusing in fluid flow through spiral ducts: dynamics, tipping phenomena and particle separation. *J. Fluid Mech.* 990, A13 (2024). https://doi.org/10.1017/jfm.2024.487
- Di Carlo D. Inertial microfluidics. *Lab Chip* 9, 3038 (2009). https://doi.org/10.1039/b912547g

## How this was built

I developed this with an AI coding assistant (Claude). I directed the scientific scope, chose the papers to validate against and reviewed each audit. The assistant wrote the code and ran the validation harness.

## Licence

MIT. See [LICENSE](LICENSE).
