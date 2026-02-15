python -m src.analysis.cli.analyze_incremental --step collect --force
python -m src.analysis.cli.analyze_incremental --step metrics --force
python -m src.analysis.cli.compute_mmgbsa_safe --force --snapshots 100 --sample-window-ns 1.0 --timestep-fs 2.0 --workers 8
python -m src.analysis.cli.analyze_incremental --step plots --force
