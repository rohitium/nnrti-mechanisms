#!/bin/bash
#
# Back up irreplaceable raw pmx NEQ FEP data from /scratch to $GROUP_HOME.
#
# WHY
#   /scratch is purged after 90 days. Raw switch data cannot be regenerated
#   without re-running the campaign, so it is copied to $GROUP_HOME, which is
#   durable.
#
#   docs/methods/md-run-data.md already mandates this: "Raw heavy data is not in git and not
#   only on scratch ... durable storage ($GROUP_HOME) with sha256 manifests."
#
# WHAT IT DOES
#   Per leg, writes one .tar.gz to $DEST, records sha256 + byte count + file
#   count in a manifest CSV, then verifies the archive is readable. Source files
#   are only ever READ -- nothing on scratch is moved or deleted.
#
# TIERS
#   TIER=1 (default)  dgdl.xvg + equil.gro + system.{gro,top} + itp + mdp +
#                     manifests + analysis/. The scientific core: every free
#                     energy in the paper can be recomputed from this. Small,
#                     because dgdl.xvg is text and compresses hard.
#   TIER=2            adds equil.trr / equil.tpr -- the endpoint trajectories.
#                     Large (this is most of the 380 GB), but it is what lets
#                     you re-extract snapshots and re-run switches WITHOUT
#                     redoing equilibration.
#   TIER=3            tier 1 + equil.xtc / .tpr / .cpt + the switch seed
#                     snapshots, but NOT equil.trr. The .trr is full precision
#                     with velocities and is ~2.0 GB per endpoint against the
#                     .xtc's 313 MB; dropping it removes ~120 GB while keeping
#                     every ensemble that is actually analysed. This is the
#                     tier chosen for the JCIM revision archive (2026-10).
#
# USAGE
#   bash ops/slurm/cluster/backup_fep_raw.sh                 # tier 1, 4 idle legs
#   TIER=2 bash ops/slurm/cluster/backup_fep_raw.sh          # + endpoint trajectories
#   LEGS="wt_to_G190E" bash ops/slurm/cluster/backup_fep_raw.sh   # after G190E finishes
#
#   This ALWAYS RUNS AS A BATCH JOB. Invoked from a login node it submits itself
#   via sbatch and returns immediately -- compression is real work and login
#   nodes are shared, throttled, and kill long processes. Set ALLOW_LOGIN_NODE=1
#   only for a deliberate small test.
#
#   Interrupted runs are safe to re-run: archives are written to .part and only
#   renamed on success, and the skip check requires a verified manifest entry.
#
# SAFETY
#   Never back up a leg with jobs actively writing to it -- tar would capture
#   half-written files. The default LEGS list deliberately EXCLUDES wt_to_G190E,
#   which is being rebuilt at 20 ns equilibration. The script refuses to run if
#   it sees queued/running pmx jobs unless FORCE=1.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"
# Walk up to the repository root rather than assuming a fixed depth: this file
# lives at ops/slurm/cluster/ in the current layout and scripts/sherlock/ in the
# older checkout still on Sherlock.
if [ -z "${PROJECT_ROOT:-}" ]; then
    PROJECT_ROOT="$SCRIPT_DIR"
    while [ "$PROJECT_ROOT" != "/" ] && [ ! -d "$PROJECT_ROOT/.git" ]; do
        PROJECT_ROOT="$(dirname "$PROJECT_ROOT")"
    done
    [ -d "$PROJECT_ROOT/.git" ] || { echo "cannot locate repository root from $SCRIPT_DIR" >&2; exit 1; }
fi
cd "$PROJECT_ROOT"

: "${GROUP_HOME:?GROUP_HOME is not set -- are you on Sherlock?}"

TIER="${TIER:-1}"
STAMP="${STAMP:-$(date +%Y-%m-%d)}"
DEST="${DEST:-$GROUP_HOME/nnrti-mechanisms/fep_raw/$STAMP}"
LEGS="${LEGS:-wt_to_K103N K103N_to_K103N_M230L K103N_to_K103N_P225H K103N_to_L100I_K103N}"
LEGS_ROOT="results/analysis/fep_pmx/legs"
MANIFEST="$DEST/BACKUP_MANIFEST.csv"

# --- self-submit; compression does not belong on a shared login node ----------
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${ALLOW_LOGIN_NODE:-0}" != "1" ]; then
    CPUS="${SBATCH_CPUS:-8}"
    if [ "$TIER" -ge 2 ]; then WALL="${SBATCH_TIME:-24:00:00}"; else WALL="${SBATCH_TIME:-04:00:00}"; fi
    mkdir -p logs/fep_backup
    echo "Submitting as a batch job (tier $TIER, ${CPUS} cpus, ${WALL} wall)."
    echo "Login nodes are shared and throttled; set ALLOW_LOGIN_NODE=1 to override."
    JOB=$(sbatch --parsable \
        -p "${SBATCH_PARTITION:-normal}" -t "$WALL" -c "$CPUS" --mem "${SBATCH_MEM:-8G}" \
        -J "fep_backup_t${TIER}" \
        -o "logs/fep_backup/fep_backup_t${TIER}_%j.log" \
        --wrap "cd '$PROJECT_ROOT' && TIER='$TIER' STAMP='$STAMP' LEGS='$LEGS' DEST='$DEST' PIGZ_THREADS='$CPUS' bash '$SELF'")
    echo "  job $JOB"
    echo "  log logs/fep_backup/fep_backup_t${TIER}_${JOB}.log"
    echo "  watch: tail -f logs/fep_backup/fep_backup_t${TIER}_${JOB}.log"
    exit 0
fi

echo "=========================================="
echo "pmx NEQ raw-data backup"
echo "=========================================="
echo "Source : $PROJECT_ROOT/$LEGS_ROOT"
echo "Dest   : $DEST"
echo "Tier   : $TIER"
echo "Legs   : $LEGS"
echo ""

# --- refuse to archive a leg that is being written to -------------------------
if [ "${FORCE:-0}" != "1" ]; then
    ACTIVE="$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -c "pmx_neq" || true)"
    if [ "${ACTIVE:-0}" -gt 0 ]; then
        for leg in $LEGS; do
            if [ "$leg" = "wt_to_G190E" ]; then
                echo "REFUSING: $ACTIVE pmx_neq job(s) are queued/running and wt_to_G190E" >&2
                echo "is in the leg list. Archiving a leg mid-write captures torn files." >&2
                echo "Wait for the campaign to finish, or set FORCE=1 if you are certain." >&2
                exit 1
            fi
        done
        echo "NOTE: $ACTIVE pmx_neq job(s) active, but none of the selected legs is"
        echo "      wt_to_G190E, so the selected legs are quiescent. Proceeding."
        echo ""
    fi
fi

mkdir -p "$DEST"

# --- progress reporting -------------------------------------------------------
# STATUS.txt answers "was it submitted / is it running / did it finish"; the
# PROGRESS.log is append-only so it can be tailed from anywhere.
STATUS_FILE="$DEST/STATUS.txt"
PROGRESS_LOG="$DEST/PROGRESS.log"
HEARTBEAT_SECONDS="${HEARTBEAT_SECONDS:-60}"

_stamp() { date -u +%Y-%m-%dT%H:%M:%SZ; }
_progress() { echo "$(_stamp) $*" >> "$PROGRESS_LOG"; }
_status() {
    printf 'state=%s\njob_id=%s\nnode=%s\ntier=%s\ndest=%s\nupdated=%s\n' \
        "$1" "${SLURM_JOB_ID:-none}" "$(hostname)" "$TIER" "$DEST" "$(_stamp)" > "$STATUS_FILE"
}
trap '_status FAILED; _progress "FAILED (exit $?)"' ERR

_status RUNNING
_progress "=== tier $TIER backup started: job ${SLURM_JOB_ID:-none} on $(hostname)"
_progress "legs: $LEGS"

if [ ! -f "$MANIFEST" ]; then
    echo "leg,tier,archive,sha256,bytes,files,source_bytes,created_utc" > "$MANIFEST"
fi

TOTAL_IN=0
TOTAL_OUT=0

for leg in $LEGS; do
    SRC="$LEGS_ROOT/$leg"
    if [ ! -d "$SRC" ]; then
        echo "SKIP $leg -- no such directory: $SRC" >&2
        continue
    fi

    OUT="$DEST/${leg}_tier${TIER}.tar.gz"
    PART="$OUT.part"
    # Skip only if the archive exists AND is recorded in the manifest. A bare
    # existence check would treat a truncated archive from an interrupted run as
    # complete -- silent data loss of exactly the kind this script exists to
    # prevent. Interrupted runs leave a .part, which is not a valid skip.
    if [ -f "$OUT" ] && grep -q ",$(basename "$OUT")," "$MANIFEST" 2>/dev/null; then
        echo "SKIP $leg -- archive exists and is verified in the manifest"
        continue
    fi
    if [ -f "$OUT" ]; then
        echo "WARN $leg -- archive exists but is NOT in the manifest (interrupted run?)."
        echo "     Re-creating it. Previous file moved to $OUT.unverified"
        mv "$OUT" "$OUT.unverified"
    fi
    rm -f "$PART"

    LIST="$(mktemp)"
    trap 'rm -f "$LIST"' EXIT

    # Tier 1: the scientific core.
    find "$SRC" \( \
        -name 'dgdl.xvg' -o \
        -name 'equil.gro' -o \
        -name 'system.gro' -o \
        -name 'system.top' -o \
        -name '*.itp' -o \
        -name '*.mdp' -o \
        -name 'neq_manifest.csv' -o \
        -name 'neq_prepare.json' -o \
        -name 'residue_map.json' -o \
        -name 'mutation.script' -o \
        -name 'analysis.json' -o \
        -name 'integ_*.dat' -o \
        -name 'results.txt' \
        \) -type f -print > "$LIST"

    case "$TIER" in
        2)
            find "$SRC" \( -name 'equil.trr' -o -name 'equil.tpr' -o -name 'equil.xtc' \) \
                -type f -print >> "$LIST"
            ;;
        3)
            # Endpoint ensembles WITHOUT the redundant full-precision .trr,
            # plus the snapshots that seeded each switch.
            find "$SRC" \( -name 'equil.xtc' -o -name 'equil.tpr' -o -name 'equil.cpt' \) \
                -type f -print >> "$LIST"
            find "$SRC" -path '*/neq/snapshots/*' -type f -print >> "$LIST"
            ;;
    esac

    NFILES="$(wc -l < "$LIST" | tr -d ' ')"
    if [ "$NFILES" -eq 0 ]; then
        echo "SKIP $leg -- no matching files" >&2
        rm -f "$LIST"; continue
    fi

    SRC_BYTES="$(xargs -a "$LIST" du -cb 2>/dev/null | tail -1 | cut -f1 || echo 0)"

    echo "-> $leg : $NFILES files, $(numfmt --to=iec "$SRC_BYTES" 2>/dev/null || echo "$SRC_BYTES B") uncompressed"
    _progress "START $leg : $NFILES files, $(numfmt --to=iec "$SRC_BYTES" 2>/dev/null || echo "$SRC_BYTES B") uncompressed"
    LEG_T0=$(date +%s)

    # pigz parallelises gzip across cores; dgdl.xvg is text, so -1 already gets
    # most of the ratio at a fraction of the CPU. Falls back to plain gzip.
    # The tar runs in the background so the heartbeat can report the growing
    # archive -- a tier-3 leg takes long enough that silence looks like a hang.
    if command -v pigz >/dev/null 2>&1; then
        ( tar -cf - -T "$LIST" | pigz -1 -p "${PIGZ_THREADS:-8}" > "$PART" ) &
    else
        ( tar -cf - -T "$LIST" | gzip -1 > "$PART" ) &
    fi
    TAR_PID=$!
    while kill -0 "$TAR_PID" 2>/dev/null; do
        sleep "$HEARTBEAT_SECONDS"
        kill -0 "$TAR_PID" 2>/dev/null || break
        NOW_BYTES="$(stat -c %s "$PART" 2>/dev/null || echo 0)"
        _progress "  .. $leg $(numfmt --to=iec "$NOW_BYTES" 2>/dev/null || echo "${NOW_BYTES}B") written after $(( $(date +%s) - LEG_T0 ))s"
    done
    wait "$TAR_PID"
    mv "$PART" "$OUT"

    SHA="$(sha256sum "$OUT" | cut -d' ' -f1)"
    OUT_BYTES="$(stat -c %s "$OUT")"

    # Verify the archive is readable and complete before trusting it.
    TAR_N="$(tar -tzf "$OUT" | wc -l | tr -d ' ')"
    if [ "$TAR_N" -ne "$NFILES" ]; then
        echo "   FAILED verification: archive holds $TAR_N entries, expected $NFILES" >&2
        mv "$OUT" "$OUT.CORRUPT"
        rm -f "$LIST"; continue
    fi

    echo "$leg,$TIER,$(basename "$OUT"),$SHA,$OUT_BYTES,$NFILES,$SRC_BYTES,$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$MANIFEST"
    _progress "DONE  $leg : $(numfmt --to=iec "$OUT_BYTES" 2>/dev/null || echo "$OUT_BYTES B") compressed, $TAR_N entries, $(( $(date +%s) - LEG_T0 ))s"
    echo "   OK  $(numfmt --to=iec "$OUT_BYTES" 2>/dev/null || echo "$OUT_BYTES B") compressed, $TAR_N entries verified"
    echo "   sha256 $SHA"

    TOTAL_IN=$((TOTAL_IN + SRC_BYTES))
    TOTAL_OUT=$((TOTAL_OUT + OUT_BYTES))
    rm -f "$LIST"
done

_status DONE
_progress "=== tier $TIER backup finished"

echo ""
echo "=========================================="
echo "Backed up : $(numfmt --to=iec "$TOTAL_IN" 2>/dev/null || echo "$TOTAL_IN B") -> $(numfmt --to=iec "$TOTAL_OUT" 2>/dev/null || echo "$TOTAL_OUT B")"
echo "Manifest  : $MANIFEST"
echo ""
echo "Copy the manifest into the repo and commit it (it is small, and it is the"
echo "record of what exists in durable storage):"
echo "  cp $MANIFEST manifests/fep_raw_backup_$STAMP.csv && git add -f manifests/fep_raw_backup_$STAMP.csv"
echo ""
echo "Verify later with:"
echo "  cd $DEST && sha256sum -c <(awk -F, 'NR>1{print \$4\"  \"\$3}' BACKUP_MANIFEST.csv)"
echo "=========================================="
