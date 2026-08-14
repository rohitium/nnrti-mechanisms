# Occupancy statistics (Atanu #5)

Replicate mean ± 95% CI (Student-t, n=3) and Welch t-test vs WT for draft-highlighted
residue contacts. Global FWER p from the existing trajectory-label permutation null.

## Verdict key

- **supported** — exceeds global 90th-percentile null **and** Welch p < 0.05
- **mixed** — passes exactly one of the two
- **descriptive_only** — fails both; do not treat as significant in main text

## Callouts for draft

### supported

- **V106A+L234I / SER105**: Δ = +0.878 [+0.806, +0.950]; Welch p = 0.000119; global FWER p = 0.0004
- **G190E / VAL179**: Δ = -0.842 [-1.009, -0.674]; Welch p = 0.000543; global FWER p = 0.0005
- **V106A+F227L / SER105**: Δ = +0.820 [+0.609, +1.032]; Welch p = 0.00328; global FWER p = 0.0047
- **V106A / SER105**: Δ = +0.751 [+0.607, +0.895]; Welch p = 0.00161; global FWER p = 0.0284
- **V106A+P225H / SER105**: Δ = +0.747 [+0.372, +1.123]; Welch p = 0.0131; global FWER p = 0.0517

### mixed

- **V106A+L234I / LYS104**: Δ = +0.672 [+0.609, +0.735]; Welch p = 0.000542; global FWER p = 0.179
- **V106A+P225H / LYS104**: Δ = +0.597 [+0.176, +1.018]; Welch p = 0.0176; global FWER p = 0.447
- **V106A+F227L / LYS104**: Δ = +0.534 [+0.317, +0.750]; Welch p = 0.00172; global FWER p = 0.58
- **V106A / LYS104**: Δ = +0.393 [+0.009, +0.776]; Welch p = 0.0353; global FWER p = 0.927

### descriptive_only

- **V106I+F227C / RES227**: Δ = -0.261 [-0.637, +0.115]; Welch p = 0.0865; global FWER p = 1
- **Y188L / LYS102**: Δ = +0.188 [-0.236, +0.611]; Welch p = 0.195; global FWER p = 1
- **G190S / VAL179**: Δ = -0.143 [-0.275, -0.012]; Welch p = 0.114; global FWER p = 1
- **Y188L / PRO225**: Δ = -0.133 [-0.902, +0.637]; Welch p = 0.551; global FWER p = 1
- **V106I+F227C / SER105**: Δ = +0.125 [-0.124, +0.374]; Welch p = 0.162; global FWER p = 1
- **K103N+P225H / PRO225**: Δ = -0.108 [-0.383, +0.167]; Welch p = 0.356; global FWER p = 1
- **Y188L / TYR188**: Δ = -0.024 [-0.481, +0.432]; Welch p = 0.889; global FWER p = 1
- **G190A / VAL179**: Δ = -0.023 [-0.255, +0.208]; Welch p = 0.785; global FWER p = 1
- **V106A+F227L / RES227**: Δ = -0.012 [-0.117, +0.093]; Welch p = 0.748; global FWER p = 1
- **V106I / SER105**: Δ = +0.008 [-0.049, +0.066]; Welch p = 0.609; global FWER p = 1

## Draft language guidance

- Keep Ser105 (V106A family) and Val179 (G190E) as load-bearing reporters if supported/mixed.
- Downgrade Y188L Lys102 / Pro225 and V106I+F227C Phe227 if descriptive_only — phrase as ‘observed shift, not significant under replicate tests’.
- n = 3 CIs are wide; never write ‘p confirms mechanism’.
