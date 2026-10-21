# Theta–gamma phase–amplitude coupling from a second-order Zwanzig–Mori equation

Code and figures for the paper

> Y. Luo, D. Hsu, S. Priyadarshini, J. Straub.
> *Theta–gamma phase–amplitude coupling from a second-order Zwanzig–Mori equation.*
> Plan for submit to *Chaos, Solitons & Fractals* (2026).

Of the ten second-order theta–gamma terms of the Zwanzig–Mori equation, only the
velocity–velocity term g<sub>vv</sub> v<sub>θ</sub> v<sub>γ</sub> locks gamma to the theta peak. The resulting model

```
x_θ'' = −ω_θ² x_θ − G_θ0 x_θ' + σ ξ_θ
x_γ'' = −ω_γ² x_γ − (G_γ0 − g_vv x_θ') x_γ' + σ ξ_γ
```

has a closed-form modulation index MI(κ), κ = g<sub>vv</sub> √(2Θ<sub>θ</sub>) / (2ω<sub>θ</sub>), and a stability boundary
g<sub>crit</sub> ≈ 2ω<sub>γ</sub> / (3√Θ<sub>θ</sub>).

## Contents

```
matlab/    core functions (MATLAB R2021a or newer, no toolboxes) and two demos
figures/   the figures of the paper (PDF and PNG)
```

| File | Purpose |
|---|---|
| `tg_sim.m` | simulator of the model (semi-implicit step with the coupling inside the trapezoidal friction step) |
| `tg_sim_quadratic.m` | gamma driven by the general quadratic theta–gamma terms |
| `tg_streams.m`, `tg_noise_chunk.m`, `tg_const.m` | random streams per seed and model constants |
| `tg_pac_record.m` | PAC of one record: Hilbert phase and envelope, log-linear detrending, Tort MI |
| `tg_surrogate_test.m` | within-cycle circular-shift surrogates |
| `tg_mi_vm.m` | closed-form MI of a von Mises envelope |
| `tg_lyapunov_grid.m` | largest Lyapunov exponent of the gamma perturbation along stochastic theta paths |
| `tg_floquet2.m`, `tg_lambda2.m`, `tg_adiabatic_gstar.m` | fixed-amplitude Floquet growth rate and its adiabatic average |
| `demo_pac.m` | MI of simulated records against the closed form |
| `demo_boundary.m` | closed-form, adiabatic and direct estimates of the stability boundary |

## Quick start

```matlab
cd matlab
demo_pac
demo_boundary
```

## Citation

See `CITATION.cff`. Please cite the paper when using this code.

## License

MIT; see `LICENSE`.
