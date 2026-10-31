# Apo MD runs — how to read these directories

The apo (ligand-free) replicates exist to answer Reviewer 1's question about
whether mutation-induced remodelling occurs in both apo and holo RT. Target is
100 ns per replicate: **50,000,000 steps at a 2 fs timestep**
(`timestep_fs = 2.0`, `production_steps_target = production_ns * 1e6 /
timestep_fs`, `src/nnrti/md/openmm/md_protocol.py`).

## Do not trust the run JSON for run length

`*_apo_repNN.json` is metadata, not a measurement.

| field | why it misleads |
|---|---|
| `status` | reads `"ok"` on truncated runs — v106a rep_02 was `ok` at 5,000,000 of 50,000,000 steps |
| `md_production_steps_completed` | wrong in *both* directions — g190e rep_02 claimed 100 ns against a 92 ns trajectory; v106i_f227c rep_03 claimed 10 ns while holding 64 ns |
| the file itself | written only on a clean exit, so a job killed by the Slurm wall clock leaves a full trajectory and **no JSON at all** |

That last point matters most: on 2026-10-31 it made 19 of 60 replicates look
empty when each held 67–95 ns.

## Authoritative sources, in order

1. **`*_md_state.csv`** — the maximum value of the `Step` column. This is the
   integrator's cumulative step count and it survives checkpoint resumes.
2. **`*_md.state.json`** (`{"step": N}`) and **`*_md.state.xml`**
   (`stepCount`, and `time` in ps). OpenMM's own record, independent of (1),
   so the two together are a genuine cross-check.
3. **Absence of `*_md.state.*`** means the run finished cleanly and its resume
   state was cleaned up. Presence means it was interrupted at least once — but
   it may still have reached target afterwards.

```bash
# authoritative length of one replicate, in ns
awk -F, '!/^#/ && $1+0>m {m=$1} END {printf "%.1f ns\n", m*2/1e6}' \
    <rep_dir>/*_md_state.csv
```

## Resume appends — so frame counts are not a clock

Resuming appends to both the DCD and the state CSV (`append_dcd` and `append=`
in `md_protocol.py`). After any resume:

- frame counts exceed the nominal 200 — up to 380 observed — with duplicated
  frames at each segment boundary;
- the `Step` column is non-monotonic; a decrease marks a resume boundary, not a
  restart from zero.

**Never derive a time axis from frame index or frame count**, and never assume
`n_frames / 200 * 100 ns`. This is the same class of error as the DCD header
bug: derive time from the step record, not from the trajectory file.

## pbcfix sidecars go stale silently

`nnrti.cli.fix_pbc_trajectories` skips any `*_analysis_pbcfix.dcd` that already
exists, so a sidecar written before a run was extended stays stale and nothing
warns you. On 2026-10-31, 23 of them did not match their raw trajectories —
several held 7 frames against 185. Before using them, compare frame counts
against the raw `*_analysis.dcd`, and pass `--force` to rebuild (back the
originals up first; they are gitignored, so git cannot recover them).
