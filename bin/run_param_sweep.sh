#!/usr/bin/env bash
# Family-parameter sweep harness for the family-profile pathway (ADR-0002 Q8, issue #6;
# extended 2026-09-03 to cover the presence-calling coverage floor, issue TBD).
#
# Runs the mmseqs family-profile pathway across the parameter grid
#   min-seq-id {0.2,0.3,0.4,0.5} x cov {0.5,0.7} x hmm-presence-E {1e-3,1e-5}
#   x hmm-presence-cov {0.5,0.3} x hmm-presence-min-residues {0,100}
# scores each grid point on #families / #novelties / BUSCO clustering recovery / BUSCO
# PRESENCE recovery / control recall+FP / TBLASTN-removal, and collates one scored table
# so the shipped default can be picked at the knee (bin/collate_sweep.py) rather than
# guessed.
#
# This is an operational launcher (like run_pezizo5.sh) — it submits real cluster jobs and
# is NOT run automatically. Edit the CONFIG block, then: bash bin/run_param_sweep.sh
#
# NOTE on grid size: the full cross product of all five dimensions below is
# 4x2x2x2x2 = 128 runs. --family_min_seq_id/--family_cov (clustering) and
# --hmm_presence_evalue/--hmm_presence_cov/--hmm_presence_min_residues (presence-calling)
# are two largely independent concerns (the latter act downstream of clustering, on
# already-built family HMMs) -- when specifically tuning the presence-calling side,
# shrink MIN_SEQ_IDS/COVS to a single element (the already-chosen clustering default)
# rather than running the full cross product.
#
# Prerequisites:
#   * BUSCO fungi_odb12 full_table.tsv per proteome, listed in BUSCO_TABLES as
#     SHORT=path (curation-free family-clustering quality, bin/busco_family_recovery.py).
#   * BUSCO_OUTGROUP_TABLES: BUSCO full_table.tsv for OUTGROUP (non-seed-group) proteomes,
#     also SHORT=path. Required for the presence-recovery metric
#     (bin/busco_presence_recovery.py) to produce a non-degenerate result -- a seed-group
#     species' own presence is trivial by clustering construction (see that script's
#     docstring), so BUSCO_TABLES alone (ingroup) cannot test hmm_presence_cov/
#     hmm_presence_min_residues at all. BUSCO_TABLES and BUSCO_OUTGROUP_TABLES may overlap
#     with BUSCO_TABLES for species that ARE seed-group members without harm (they simply
#     contribute no scorable pairs), but at least one genuinely non-member species is
#     needed for a real signal.
#   * A curated controls CSV (configs/controls/<clade>.controls.csv) for recall/FP
#     (optional — recall/FP columns are left blank if omitted).
set -euo pipefail

# ------------------------------------------------------------------ CONFIG (edit me)
CONFIG=${CONFIG:-configs/pezizo5.csv}
DATA_DIR=${DATA_DIR:-data}
CLADE=${CLADE:-pezizo5}
PFAM=${PFAM:-db/pfam/38.2/Pfam-A.hmm}
SWISSPROT=${SWISSPROT:-db/uniprot/uniprot_sprot.fasta.dmnd}
MODELORGS=${MODELORGS:-configs/modelorgs.yaml}
CONTROLS=${CONTROLS:-configs/controls/${CLADE}.controls.csv}
# NXF_SITE_CONFIG defaults to conf/ucr_hpcc_slurm.config (2026-09-04 fix): without it,
# `-profile slurm` alone gets none of that file's queue routing, AVX2 node constraints,
# retry/memory scaling, or BUILD_CHUNK preempt-queue placement -- nextflow.config's own
# minimal `slurm{}` profile block just sets executor=slurm with no queue at all, so SLURM's
# own unspecified-partition default applies. Confirmed causing real failures at deep_broad_1kfg
# scale: a BUILD_CHUNK task killed by preemption, and an HMMSEARCH_CHUNK task OOM-killed on
# all 3 attempts (a flat memory request with no site-specific per-attempt scaling). Set to
# empty to opt out (e.g. running off-cluster).
PROFILE=${PROFILE:--profile slurm}
NXF_SITE_CONFIG=${NXF_SITE_CONFIG:-conf/ucr_hpcc_slurm.config}
OUTDIR=${OUTDIR:-results}
SWEEP_DIR=${SWEEP_DIR:-results/sweep_${CLADE}}
# BUSCO full_table.tsv per proteome, space-separated SHORT=path (edit for your run):
BUSCO_TABLES=${BUSCO_TABLES:-}
# BUSCO full_table.tsv for OUTGROUP proteomes -- see "Prerequisites" above. Required for
# a non-degenerate presence_recovery number.
BUSCO_OUTGROUP_TABLES=${BUSCO_OUTGROUP_TABLES:-}
# Any FAMILY_HMMSEARCH domtblout from the run, for the eligible-BUSCO-vs-whole-proteome
# length-representativeness check (bin/busco_presence_recovery.py --reference-domtblout).
# Set per grid point below if left empty (first ingroup proteome's domtblout).
REFERENCE_DOMTBLOUT=${REFERENCE_DOMTBLOUT:-}

# Sweep grid (ADR-0002 Q8 + 2026-09-03 coverage-floor extension).
# Each dimension can be scoped down via a space-separated env var (e.g.
# MIN_SEQ_IDS_LIST="0.3" to hold clustering identity fixed while sweeping only the
# presence-calling dimensions below) -- see the grid-size NOTE above.
read -ra MIN_SEQ_IDS <<< "${MIN_SEQ_IDS_LIST:-0.2 0.3 0.4 0.5}"
read -ra COVS <<< "${COVS_LIST:-0.5 0.7}"
read -ra HMM_EVALUES <<< "${HMM_EVALUES_LIST:-1e-3 1e-5}"
read -ra HMM_COVS <<< "${HMM_COVS_LIST:-0.5 0.3}"
read -ra MIN_RESIDUES <<< "${MIN_RESIDUES_LIST:-0 100}"

# ---------------------------------------------------------------------------- setup
REPO_ROOT=$(pwd)
BIN=$(cd "$(dirname "$0")" && pwd)

# Isolate this sweep's `nextflow run` launch directory (2026-09-03 fix): Nextflow's
# session lock/resume-history (.nextflow/cache/<session>/db/LOCK, .nextflow/history) is
# keyed to the LAUNCH (current working) directory, not to --project/--outdir. Running two
# sweeps concurrently from the same directory made a plain `-resume` (no explicit run
# name) in one invocation try to attach to the OTHER invocation's still-running session,
# which fails outright ("Unable to acquire lock on session..."; confirmed when the
# pezizo5/sordariales_shallow/deep_broad_1kfg sweeps were launched together -- deep_broad_1kfg
# lost its entire first attempt to this). Fix: give every sweep (keyed by --clade) its own
# launch directory, and resolve every path-bearing variable to absolute BEFORE cd'ing there
# (nextflow.config alongside main.nf is still auto-loaded regardless of launch directory --
# only the launch-dir-relative .nextflow/ metadata needs isolating). DSL2 `include` paths in
# main.nf resolve relative to main.nf's own location, not the launch directory, so this is
# safe.
abspath() { python3 -c "import os,sys; print(os.path.abspath(sys.argv[1]))" "$1"; }
abspath_tables() {  # abspath every "SHORT=path" token in a space-separated string
    local out="" tok short path
    for tok in $1; do
        short="${tok%%=*}"; path="${tok#*=}"
        out="$out ${short}=$(abspath "$path")"
    done
    echo "$out"
}
CONFIG=$(abspath "$CONFIG")
DATA_DIR=$(abspath "$DATA_DIR")
PFAM=$(abspath "$PFAM")
SWISSPROT=$(abspath "$SWISSPROT")
MODELORGS=$(abspath "$MODELORGS")
CONTROLS=$(abspath "$CONTROLS")
[ -n "$NXF_SITE_CONFIG" ] && NXF_SITE_CONFIG=$(abspath "$NXF_SITE_CONFIG")
OUTDIR=$(abspath "$OUTDIR")
mkdir -p "$SWEEP_DIR"
SWEEP_DIR=$(abspath "$SWEEP_DIR")
BUSCO_TABLES=$(abspath_tables "$BUSCO_TABLES")
BUSCO_OUTGROUP_TABLES=$(abspath_tables "$BUSCO_OUTGROUP_TABLES")
[ -n "$REFERENCE_DOMTBLOUT" ] && REFERENCE_DOMTBLOUT=$(abspath "$REFERENCE_DOMTBLOUT")

LAUNCH_DIR="$REPO_ROOT/.nf_launch/${CLADE}"
mkdir -p "$LAUNCH_DIR"
cd "$LAUNCH_DIR"

METRICS="$SWEEP_DIR/sweep_metrics.tsv"
printf 'min_seq_id\tcov\thmm_evalue\thmm_cov\thmm_residues\tn_families\tn_novelties\tbusco_recovery\tpresence_recovery\trecall\tfp_rate\ttblastn_removed\trun_ok\n' > "$METRICS"

count_lines() { [ -s "$1" ] && grep -cve '^\s*$' "$1" || echo 0; }

# --------------------------------------------------------------------------- sweep
for id in "${MIN_SEQ_IDS[@]}"; do
  for cov in "${COVS[@]}"; do
    for ev in "${HMM_EVALUES[@]}"; do
      for hcov in "${HMM_COVS[@]}"; do
        for res in "${MIN_RESIDUES[@]}"; do
          tag="id${id}_c${cov}_e${ev}_hc${hcov}_r${res}"
          project="sweep_${CLADE}_${tag}"
          rundir="$OUTDIR/$project"
          echo ">>> grid point $tag  (project $project)"

          # Annotation (Pfam/SwissProt) is skipped by default (ANNOTATE_SWEEP=0): none of
          # the sweep metrics below need it (matrix falls back to the un-annotated
          # presence_matrix.tsv), it is the most expensive step, and its MPI hmmsearch has
          # failed with SLURM slot-allocation errors independent of any parameter under
          # test here. Set ANNOTATE_SWEEP=1 to include it (e.g. to also inspect annotated
          # output by eye), at the cost of exposing the sweep to that failure mode.
          annotate_args=()
          if [ "${ANNOTATE_SWEEP:-0}" = "1" ]; then
            annotate_args=(--pfam_hmm "$PFAM" --swissprot_dmnd "$SWISSPROT" --modelorgs_config "$MODELORGS")
          fi

          # A single grid point's pipeline failure must not abort the whole sweep (the
          # script runs under `set -e`) -- record it (blank/worst-value metrics) and move
          # on to the next point instead.
          run_ok=1
          site_config_args=()
          [ -n "$NXF_SITE_CONFIG" ] && site_config_args=(-c "$NXF_SITE_CONFIG")
          nextflow run "$REPO_ROOT/main.nf" -resume $PROFILE "${site_config_args[@]}" \
            --config "$CONFIG" --data_dir "$DATA_DIR" \
            --run_tool phmmer --cluster_tool mmseqs \
            --family_min_seq_id "$id" --family_cov "$cov" --hmm_presence_evalue "$ev" \
            --hmm_presence_cov "$hcov" --hmm_presence_min_residues "$res" \
            "${annotate_args[@]}" \
            --outdir "$OUTDIR" --project "$project" || run_ok=0
          if [ "$run_ok" = "0" ]; then
            echo ">>> grid point $tag FAILED -- recording partial/blank metrics and continuing" >&2
          fi

          cluster_tsv="$rundir/families/families_cluster.tsv"
          families_tsv="$rundir/families/families.tsv"
          matrix="$rundir/presence_matrix.function.tsv"
          [ -s "$matrix" ] || matrix="$rundir/presence_matrix.tsv"

          n_families=$(($(count_lines "$families_tsv") > 0 ? $(count_lines "$families_tsv") - 1 : 0))
          n_novelties=$(count_lines "$rundir/candidates.txt")

          # BUSCO single-copy family-clustering recovery (curation-free, ADR-0002 Q8).
          # Insensitive to hmm_cov/hmm_residues -- see bin/busco_family_recovery.py.
          busco_recovery=""
          if [ -n "$BUSCO_TABLES" ] && [ -s "$families_tsv" ]; then
            "$BIN/busco_family_recovery.py" --tables $BUSCO_TABLES \
              --cluster-tsv "$cluster_tsv" --families "$families_tsv" \
              --output "$SWEEP_DIR/busco_${tag}.tsv" || true
            busco_recovery=$(awk -F'\t' '$1=="recovery_rate"{print $2}' \
              "$SWEEP_DIR/busco_${tag}.summary.tsv" 2>/dev/null | tr -d '\r' || true)
          fi

          # BUSCO PRESENCE recovery (curation-free, 2026-09-03): the metric that is
          # actually sensitive to hmm_cov/hmm_residues -- see bin/busco_presence_recovery.py
          # and its docstring's note that BUSCO_TABLES alone cannot test this (needs
          # OUTGROUP coverage).
          presence_recovery=""
          if [ -n "$BUSCO_TABLES" ] && [ -n "$BUSCO_OUTGROUP_TABLES" ] && [ -s "$matrix" ]; then
            ref="$REFERENCE_DOMTBLOUT"
            if [ -z "$ref" ]; then
              ref=$(ls "$rundir"/family_hmmsearch/*.domtblout 2>/dev/null | head -1 || true)
            fi
            ref_arg=""
            [ -n "$ref" ] && [ -s "$ref" ] && ref_arg="--reference-domtblout $ref"
            "$BIN/busco_presence_recovery.py" --tables $BUSCO_TABLES $BUSCO_OUTGROUP_TABLES \
              --cluster-tsv "$cluster_tsv" --families "$families_tsv" --matrix "$matrix" \
              $ref_arg \
              --output "$SWEEP_DIR/presence_${tag}.tsv" || true
            presence_recovery=$(awk -F'\t' '$1=="presence_recovery_rate"{print $2}' \
              "$SWEEP_DIR/presence_${tag}.summary.tsv" 2>/dev/null | tr -d '\r' || true)
          fi

          # Control recall / false-novelty rate.
          recall=""; fp_rate=""
          if [ -s "$CONTROLS" ] && [ -s "$matrix" ]; then
            "$BIN/score_controls.py" --controls "$CONTROLS" --matrix "$matrix" \
              --cluster-tsv "$cluster_tsv" --families "$families_tsv" --config "$CONFIG" \
              --profiles "$rundir/families/family_profiles.hmm" \
              --output "$SWEEP_DIR/controls_${tag}.tsv" || true
            recall=$(awk -F'\t' '$1=="recall"{print $2}' \
              "$SWEEP_DIR/controls_${tag}.summary.tsv" 2>/dev/null | tr -d '\r' || true)
            fp_rate=$(awk -F'\t' '$1=="fp_rate"{print $2}' \
              "$SWEEP_DIR/controls_${tag}.summary.tsv" 2>/dev/null | tr -d '\r' || true)
          fi

          # TBLASTN-removal count (annotation-artifact rate): candidates with any outgroup
          # genomic hit — reported as a quality signal, not a filter here.
          tblastn_removed=0
          tbsum="$rundir/tblastn_summary.tsv"
          if [ -s "$tbsum" ]; then
            tblastn_removed=$(awk -F'\t' 'NR>1{for(i=2;i<=NF;i++) if($i==1){c++;break}} END{print c+0}' "$tbsum")
          fi

          printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$id" "$cov" "$ev" "$hcov" "$res" "$n_families" "$n_novelties" \
            "$busco_recovery" "$presence_recovery" "$recall" "$fp_rate" "$tblastn_removed" \
            "$run_ok" >> "$METRICS"
        done
      done
    done
  done
done

# --------------------------------------------------------------------------- collate
"$BIN/collate_sweep.py" --metrics "$METRICS" --output "$SWEEP_DIR/sweep_scores.tsv"
echo ">>> scored table: $SWEEP_DIR/sweep_scores.tsv"
echo ">>> wire the recommended default into nextflow.config (family_min_seq_id / family_cov / hmm_presence_evalue / hmm_presence_cov / hmm_presence_min_residues)"
