# Snakemake Prototype

This is a minimal Snakemake scaffold that wraps the existing scripts without changing their behavior.
It is meant for development on the desktop and execution on the cluster.

## Quick start (dry run)

- `snakemake -n --use-conda`

## Notes

- Update paths in `config/config.yaml` to the cluster locations before running.
- Each rule currently produces a marker file under `results/markers/` after the script finishes.
- The scripts currently hardcode many paths; those can be parameterized later and wired through Snakemake.
- Singularity can be enabled later with `--use-singularity` once containers are defined.
