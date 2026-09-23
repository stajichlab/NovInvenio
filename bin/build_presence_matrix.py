#!/usr/bin/env python3
"""Build a protein × proteome presence/absence matrix and emit candidate IDs.

Presence scoring uses one flat significance filter plus one paralog-aware filter,
with an optional absolute-evalue override on the second:

  1. Significance filter: hit e-value < --default-evalue (flat, 1e-5).

  2. Paralog-competition filter: if the query's paralog hits the same target
     with a *better* (lower) e-value than the query does, the hit is
     disqualified.  The logic is that such a hit is better explained by the
     conserved domain shared with the paralog than by the query protein itself.
     --paralog-competition-scope controls what "the same target" means:
       proteome (default) — the paralog's best hit anywhere in the target
         *proteome* out-scores the query. Strict; can drop a real ortholog when
         the query's paralog also has a strong (but distinct) ortholog in the
         target genome — e.g. HEX-1, a derivative of eIF5A, is dropped because
         the genome's true eIF5A hits slightly harder than the HEX-1 ortholog.
       target — the paralog out-scores the query on the *same target protein*.
         Preserves a hit when the paralog wins only on a different gene, so the
         HEX-1 ortholog call survives.

     --paralog-rescue-evalue (ON by default, 1e-20 -- issue #128) puts a floor under
     filter 2: a hit is never disqualified if the query's own e-value against
     that target is already <= this threshold, regardless of how much harder the
     paralog hits. Filter 2 alone only compares the query against its paralog, so
     it can discard a hit that is independently strong evidence of homology --
     e.g. NCU11312/A7UWR3_NEUCR (a CorA-family Mg2+ transporter) hits every
     outgroup proteome at ~1e-93..1e-97, but its in-genome paralog Q7SEN2 wins
     head-to-head on every one of those same targets (~1e-98..1e-131), so filter
     2 alone zeroes out all six outgroup cells and calls it a lineage-specific
     novelty. That differs from the HEX-1/eIF5A case filter 2 is designed
     around: HEX-1 has no raw hit at all in five of six outgroups, and only a
     marginal one (~1e-9..1e-12) in the sixth -- filter 2 has essentially
     nothing to override there, so its "novel" call is a real absence-of-signal
     result, not a competition artifact. --paralog-rescue-evalue distinguishes
     the two: it keeps a hit like A7UWR3's -- strong on its own merits --
     counted as present even though its paralog wins the head-to-head, while
     leaving marginal hits like HEX-1's still subject to filter 2. Pass 0 to
     disable (no e-value is <= 0, so the rescue never fires), restoring the
     behaviour tagged baseline/pre-paralog-rescue-2026-09-20.

     --paralog-rescue-delta (optional, OFF by default) is a second rescue arm,
     OR-ed with the floor: keep a hit the paralog beat by fewer than DELTA
     orders of magnitude, i.e. log10(query_ev) - log10(paralog_ev) < DELTA. It
     asks whether the paralog explained the hit away rather than whether the hit
     is strong in absolute terms, and is scale-free where the floor is not.

     It is off by default because it measured NON-SELECTIVE on real data, stated
     plainly here so it is not re-adopted on intuition. The rescue fires per
     (protein, proteome) cell and any one rescued cell ends a novelty call, so
     the value that triggers the rule is the MINIMUM delta over a candidate's
     suppressed cells -- median 8.6 on pezizo_set1. `delta < 45` therefore fires
     for 291 of the 334 filter-2-suppressed candidates, barely narrower than
     removing all 334. Scored against TBLASTN outgroup-genome breadth >= 4 as an
     independent proxy for "the gene really is there", it removes 291 candidates
     at a 47.1% hit-rate versus 48.2% for removing every suppressed candidate
     indiscriminately -- worse than no rule -- while the 1e-20 floor removes 213
     at 57.3%. OR-ing delta onto the floor drags 57.3% down to 49.0%. Both
     controls (A7UWR3 min-delta 2.0, HEX-1 min-delta 56.4) pass under either
     arm, so they cannot distinguish them. Retained as a flag so it can be
     re-swept once alignment coverage exists (issue #129), not recommended.

     That proxy has its own unmeasured false-positive rate, so the 57.3%/48.2%
     margin is softer than it looks; see issue #135.

     Both arms are mirrored, with the same defaults and the same 0-disables
     convention, in lib/singleton_presence.py's score_singleton_hits() (the
     --cluster_tool novelty_discovery singleton branch) and in
     bin/context_presence.py's inline copy (CONTEXT_SEARCH's report-only
     NEAR_INGROUP/BROAD_OUTGROUP evidence) -- issue #138. All three carry the
     same filter 2, so a change here must be made in all three or the same
     protein gets different answers from different pathways.

  3. Other-group coverage floor (OPT-IN, OFF by default -- issue #158):
     --other-coverage-floor-qcov Q treats a hit to an *other-group* proteome (the
     absence side: outgroup in the novelty direction, ingroup in the loss
     direction) as absent when its query coverage qcov < Q percent. A narrow hit
     like spa-18's (qcov 15.1 against an unrelated 1117aa protein) is the shape of
     a domain-driven or short-motif artifact, not orthology. Identity is NOT a
     co-condition: an earlier `qcov<30 AND pident>=40` rule missed spa-18
     (pident 27.8) and rejected more true orthologs.

     Measured on pezizo_set1 (docs/superpowers/specs/
     2026-09-22-coverage-floor-sensitivity-design-handoff.md, Sec 4-5): Q=20 wrongly
     rejects 1.18% of 2367 universally single-copy BUSCO ortholog pairs and flags
     11.0% of 1231 current outgroup presence hits; Q=15 gives 0.38% and 7.7%. BUSCO
     genes are conserved core genes, so these false-rejection rates are a lower
     bound for real candidates, not an estimate. Q=15 is the chosen value when the
     floor is enabled (2026-09-23): it favours fewer false rejections.

     Query-group cells are never filtered: there a wrong rejection drops a
     candidate below --ingroup-min-frac with nothing downstream to recover it.
     The floor runs after filter 2 and its rejections are counted separately
     (stderr, and --output-coverage-floor-rejections), so a hit such as HEX-1's,
     already explained by filter 2, stays attributable to filter 2 alone.

     The floor needs the wide hit layout (issue #129: diamond/blastp qcovhsp). A
     requested floor with a missing or blank qcov on any hit it would judge (an old
     4-column cache, or any phmmer --tblout hit) is a hard error, never a silent
     no-op -- see lib/hits.py's Hit docstring.

  (2026-09-03: filter 1 used to be a per-query "paralog-cutoff" -- hit e-value
  must beat the query's own within-proteome paralog e-value, falling back to
  --default-evalue when no paralog was detected. Dropped: deriving an absolute
  significance ceiling from a self-vs-self search (reported down to E=100, not
  the usual significance floor) swung both too tight -- any protein with a
  close in-genome paralog got an unreachable cutoff, rejecting real orthologs
  regardless of relevance -- and too loose -- a protein with no real paralog
  still picked up self-search noise as its "cutoff", looser than the intended
  default. Filter 2 already does the actual paralogy test, as a direct
  head-to-head comparison instead of an absolute-magnitude proxy; filter 1 was
  redundant on top of it, not protective. See lib/singleton_presence.py's
  module docstring for the mirrored fix and empirical evidence.)

A candidate protein (from a --query-group proteome, default IN) must be:
  - present in >= --ingroup-min-frac of all --query-group proteomes
  - present in <= --other-max-frac of the other group's proteomes (default 0.0,
    i.e. absent from every other-group proteome)

--query-group OUT runs the same logic in the opposite direction — outgroup
proteomes as query, looking for genes conserved in the outgroup but absent
from the ingroup (candidate lineage-specific losses). --ingroup-min-frac is
reused as the query group's own presence threshold in both directions (pass
--ingroup-min-frac params.outgroup_min_frac from the loss-search workflow).
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / 'lib'))
from config_parser import parse_config

DEFAULT_EVALUE = 1e-5
# Filter-2 rescue floor, ON by default (issue #128). Mirrors nextflow.config's
# params.paralog_rescue_evalue. 0 disables.
DEFAULT_RESCUE_EVALUE = 1e-20

REJECTION_COLUMNS = ['query_proteome', 'query_id', 'target_proteome', 'target_id',
                     'evalue', 'qcov', 'pident', 'length', 'qlen', 'slen']


def load_paralog_info(cutoff_files):
    """Return paralog_of: protein_id -> paralog_protein_id, from all paralog_cutoffs.tsv
    files -- feeds filter 2 (paralog-competition) only; see module docstring for why the
    per-query e-value column is no longer used as a significance cutoff (filter 1)."""
    paralog_of = {}
    for path in cutoff_files:
        df = pd.read_csv(path, sep='\t')
        if df.empty:
            continue
        paralog_of.update(dict(zip(df['protein_ID'], df['paralog_protein_ID'])))
    return paralog_of


HIT_COLUMNS = ['query_id', 'target_id', 'evalue', 'bitscore',
               'query_proteome', 'target_proteome',
               # Metric columns (issue #129), carried through by bin/parse_hits.py when
               # the raw tool output has them; blank/absent otherwise. qcov feeds the
               # opt-in other-group coverage floor (filter 3, issue #158).
               'length', 'pident', 'qcov', 'scov', 'qlen', 'slen']


def load_hits(hit_files):
    """Concatenate all parsed pairwise-hit TSVs into a single DataFrame."""
    frames = []
    for tsv_path in hit_files:
        df = pd.read_csv(tsv_path, sep='\t')
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=HIT_COLUMNS)
    hits = pd.concat(frames, ignore_index=True)
    hits['evalue'] = hits['evalue'].astype(float)
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hits', nargs='+', required=True,
                    help='Parsed hit TSV files (one per query-target pair)')
    ap.add_argument('--paralog-cutoffs', nargs='+', default=[],
                    dest='paralog_cutoffs',
                    help='Per-species paralog_cutoffs.tsv files from self-vs-self search')
    ap.add_argument('--config',  required=True, help='Analysis description CSV')
    ap.add_argument('--ingroup-min-frac', type=float, default=0.75,
                    dest='ingroup_min_frac',
                    help='Presence threshold within --query-group (fraction of that '
                         'group\'s proteomes that must contain a hit)')
    ap.add_argument('--query-group', choices=['IN', 'OUT'], default='IN',
                    dest='query_group',
                    help='Which config group supplies the query proteomes candidates are '
                         'sourced from (default: IN, the novelty-search direction). OUT '
                         'runs the loss-search direction: outgroup query, ingroup must be 0.')
    ap.add_argument('--other-max-frac', type=float, default=0.0,
                    dest='other_max_frac',
                    help='Max fraction of the *other* group a candidate may still be '
                         'present in (default: 0.0 = must be absent from every other-group '
                         'proteome). Relax it in the loss direction to allow candidates '
                         'that survive in a small fraction of the ingroup — genes nearly, '
                         'but not entirely, lost.')
    ap.add_argument('--paralog-competition-scope', choices=['proteome', 'target'],
                    default='target', dest='paralog_competition_scope',
                    help="Granularity of filter 2 (paralog-competition). 'target' "
                         "(default, matches nextflow.config's pipeline default): disqualify "
                         "only if the paralog beats the query on the *same target protein* — "
                         "keeps calls where the paralog wins on a different gene (e.g. HEX-1 "
                         "vs its ancestral eIF5A). 'proteome' (original behaviour, stricter): "
                         "disqualify a hit if the query's paralog out-scores the query "
                         "anywhere in the same target proteome.")
    ap.add_argument('--default-evalue', type=float, default=DEFAULT_EVALUE,
                    dest='default_evalue',
                    help='Flat e-value significance cutoff (filter 1), applied to every hit')
    ap.add_argument('--paralog-rescue-evalue', type=float, default=DEFAULT_RESCUE_EVALUE,
                    dest='paralog_rescue_evalue',
                    help='Floor under filter 2 (paralog-competition): a hit is never '
                         'disqualified if the query\'s own e-value against that target is '
                         'already <= this threshold, no matter how much harder the paralog '
                         'hits (see module docstring for the A7UWR3-vs-HEX1 motivating case). '
                         f'Default {DEFAULT_RESCUE_EVALUE:g}, matching nextflow.config. Pass 0 '
                         'to disable and restore the pre-2026-09-20 behaviour -- no e-value is '
                         '<= 0, so the rescue never fires.')
    ap.add_argument('--paralog-rescue-delta', type=float, default=None,
                    dest='paralog_rescue_delta',
                    help='Second, OPT-IN rescue arm under filter 2, OR-ed with '
                         '--paralog-rescue-evalue: a hit is never disqualified if the paralog '
                         'beat it by fewer than this many orders of magnitude, i.e. '
                         'log10(query_evalue) - log10(paralog_evalue) < DELTA. Asks whether the '
                         'paralog explained the hit away rather than whether the hit is strong '
                         'in absolute terms. Disabled (None) by default and NOT recommended -- '
                         'measured non-selective on real data, see the module docstring.')
    ap.add_argument('--other-coverage-floor-qcov', type=float, default=None,
                    dest='other_coverage_floor_qcov',
                    help='OPT-IN filter 3: treat a hit to an other-group proteome as absent '
                         'when its query coverage (qcov, percent) is below this value. Never '
                         'applied to query-group cells. Off by default; 0 also disables. '
                         'Requires wide hit files with qcov (diamond/blastp, issue #129) -- '
                         'a missing or blank qcov is a hard error. See the module docstring '
                         'for measured rates at 15 and 20.')
    ap.add_argument('--output-coverage-floor-rejections', default=None,
                    dest='output_coverage_floor_rejections',
                    help='Optional TSV listing every hit filter 3 rejected (header only when '
                         'the floor is off). Rejections by filter 2 are never listed here.')
    ap.add_argument('--output-matrix',     required=True)
    ap.add_argument('--output-candidates', required=True)
    ap.add_argument('--output-evalues', default=None, dest='output_evalues',
                    help='Optional sidecar TSV, same shape as --output-matrix, holding the '
                         'best qualifying hit e-value per (protein, proteome) cell instead '
                         'of 0/1 (empty when absent or self-sourced) — report-only evidence, '
                         'does not affect candidate calling.')
    ap.add_argument('--output-targets', default=None, dest='output_targets',
                    help='Optional sidecar TSV, same shape as --output-matrix/--output-evalues, '
                         "holding the best qualifying hit's target protein ID per (protein, "
                         'proteome) cell (empty when absent or self-sourced) — report-only '
                         "evidence (resolved to a name via bin/extract_protein_descriptions.py's "
                         'output downstream), does not affect candidate calling.')
    args = ap.parse_args()

    samples      = parse_config(args.config)
    # Deliberately strict IN/OUT (not INGROUP_ROLES/OUTGROUP_ROLES): this is the
    # --cluster_tool pairwise matrix builder, and main.nf's SEARCH/LOSS_SEARCH
    # channels only ever feed it IN/OUT proteomes -- broadening this to the coarse
    # bands (as bin/profile_to_matrix.py's mmseqs path does) would create
    # NEAR_INGROUP/BROAD_OUTGROUP matrix columns that are always empty (never
    # actually searched), which then silently shadows the real evidence issue #48's
    # CONTEXT_SEARCH produces for those same proteomes. NEAR_INGROUP/BROAD_OUTGROUP
    # rows in a pairwise config are handled by CONTEXT_SEARCH instead (report-only,
    # candidates-only), not by this matrix.
    ingroup_ids  = {s.short for s in samples if s.group == 'IN'}
    outgroup_ids = {s.short for s in samples if s.group == 'OUT'}
    all_ids      = ingroup_ids | outgroup_ids
    query_ids    = ingroup_ids if args.query_group == 'IN' else outgroup_ids
    other_ids    = outgroup_ids if args.query_group == 'IN' else ingroup_ids

    paralog_of = load_paralog_info(args.paralog_cutoffs)

    hits = load_hits(args.hits)

    # Best evalue lookup for filter 2, keyed by how well a paralog hits the target.
    # 'proteome' scope keys on the whole target proteome; 'target' scope keys on the
    # individual target protein so a paralog only disqualifies a hit it beats head-to-head.
    if hits.empty:
        best_ev = {}
    elif args.paralog_competition_scope == 'proteome':
        best_ev = (hits.groupby(['query_proteome', 'query_id', 'target_proteome'])
                       ['evalue'].min().to_dict())
    else:  # 'target'
        best_ev = (hits.groupby(['query_proteome', 'query_id', 'target_id'])
                       ['evalue'].min().to_dict())

    # Restrict to --query-group queries, then apply both paralog-aware filters vectorised.
    ing = hits[hits['query_proteome'].isin(query_ids)].copy()

    if not ing.empty:
        # Filter 1: flat significance floor (no longer per-query paralog-derived).
        ing = ing[ing['evalue'] < args.default_evalue]

    if not ing.empty:
        # Filter 2: disqualify if the query's paralog beats it. The key is the target
        # proteome ('proteome' scope) or the individual target protein ('target' scope).
        paralog_ids = ing['query_id'].map(paralog_of)
        target_key = (ing['target_proteome'] if args.paralog_competition_scope == 'proteome'
                      else ing['target_id'])
        paralog_ev = [
            best_ev.get((qp, pid, tk)) if pid is not None and not pd.isna(pid) else None
            for qp, pid, tk in zip(ing['query_proteome'], paralog_ids, target_key)
        ]
        paralog_ev = pd.Series(paralog_ev, index=ing.index, dtype='float64')
        disqualified = paralog_ev.notna() & (paralog_ev < ing['evalue'])

        # Rescue arms, OR-ed: a hit survives filter 2 if EITHER arm covers it. They are
        # deliberately independent -- the floor asks "is this hit strong on its own?",
        # the delta arm asks "did the paralog actually explain it away?" -- so a hit
        # covered by one but not the other must still be kept. See the module docstring
        # for why only the floor is enabled by default.
        rescued = pd.Series(False, index=ing.index)
        if args.paralog_rescue_evalue:            # 0 (or None) disables this arm
            rescued |= ing['evalue'] <= args.paralog_rescue_evalue
        if args.paralog_rescue_delta is not None:
            # log10(0) is -inf, which is the correct reading at both ends: a zero-e-value
            # paralog (diamond's report for an overwhelming hit) beat the query by an
            # unbounded margin -> delta +inf -> never rescued.
            with np.errstate(divide='ignore', invalid='ignore'):
                delta = np.log10(ing['evalue']) - np.log10(paralog_ev)
            rescued |= pd.Series(delta, index=ing.index) < args.paralog_rescue_delta
        disqualified &= ~rescued
        n_filter2 = int(disqualified.sum())
        ing = ing[~disqualified]
    else:
        n_filter2 = 0

    # Filter 3 (opt-in): other-group coverage floor. Runs after filter 2 so each
    # rejected hit is attributable to exactly one mechanism.
    floor = args.other_coverage_floor_qcov
    rejected = pd.DataFrame(columns=REJECTION_COLUMNS)
    n_cells_flipped = 0
    if floor and not ing.empty:
        judged = ing['target_proteome'].isin(other_ids)
        qcov = (pd.to_numeric(ing['qcov'], errors='coerce') if 'qcov' in ing.columns
                else pd.Series(np.nan, index=ing.index))
        n_unmeasured = int((judged & qcov.isna()).sum())
        if n_unmeasured:
            sys.exit(
                f'ERROR: --other-coverage-floor-qcov {floor:g} was requested, but '
                f'{n_unmeasured} of {int(judged.sum())} other-group hits have no qcov '
                '(an old 4-column cached hit file, or a phmmer --tblout hit). Re-run the '
                'search with diamond/blast wide output (issue #129), or drop the floor. '
                'Refusing to skip those hits silently.')
        floor_hit = judged & (qcov < floor)
        rejected = ing.loc[floor_hit].reindex(columns=REJECTION_COLUMNS)
        cell_cols = ['query_proteome', 'query_id', 'target_proteome']
        cells_before = set(map(tuple, ing.loc[judged, cell_cols].to_numpy()))
        ing = ing[~floor_hit]
        cells_after = set(map(tuple, ing.loc[ing['target_proteome'].isin(other_ids),
                                             cell_cols].to_numpy()))
        n_cells_flipped = len(cells_before - cells_after)

    print(f'Hits rejected by filter 2 (paralog competition): {n_filter2} hit(s)',
          file=sys.stderr)
    if floor:
        print(f'Hits rejected by filter 3, coverage floor (qcov < {floor:g}): '
              f'{len(rejected)} hit(s); {n_cells_flipped} other-group presence cell(s) '
              'now absent', file=sys.stderr)
    if args.output_coverage_floor_rejections:
        rejected.to_csv(args.output_coverage_floor_rejections, sep='\t', index=False)

    # protein_key (query_proteome, protein_id) -> set of target proteomes with qualifying hits
    presence: dict[tuple, set] = defaultdict(set)
    for qp, qid, tp in zip(ing['query_proteome'], ing['query_id'], ing['target_proteome']):
        presence[(qp, qid)].add(tp)

    # Best (lowest) qualifying-hit e-value (and that hit's target_id) per (query_proteome,
    # query_id, target_proteome) — report-only evidence for the e-value/targets sidecars;
    # does not affect candidate calling. idxmin (not a plain groupby().min()) so hit_target
    # names the actual winning row's target, not just an independently-computed min value.
    hit_evalue: dict[tuple, float] = {}
    hit_target: dict[tuple, str] = {}
    if not ing.empty:
        best_idx = (ing.groupby(['query_proteome', 'query_id', 'target_proteome'])
                       ['evalue'].idxmin())
        best_rows = ing.loc[best_idx]
        for qp, pid, tp, ev, tid in zip(
            best_rows['query_proteome'], best_rows['query_id'], best_rows['target_proteome'],
            best_rows['evalue'], best_rows['target_id'],
        ):
            hit_evalue[(qp, pid, tp)] = ev
            hit_target[(qp, pid, tp)] = tid

    # Build the full matrix (always emit the id columns + one column per proteome,
    # so an empty result still writes a well-formed header).
    sorted_ids = sorted(all_ids)
    columns = ['protein_id', 'source_proteome'] + sorted_ids
    rows = []
    evalue_rows = []
    target_rows = []
    for (qp, pid), hit_proteomes in presence.items():
        all_present = hit_proteomes | {qp}
        row = {'protein_id': pid, 'source_proteome': qp}
        ev_row = {'protein_id': pid, 'source_proteome': qp}
        tgt_row = {'protein_id': pid, 'source_proteome': qp}
        for sp in sorted_ids:
            row[sp] = int(sp in all_present)
            ev_row[sp] = '' if sp == qp else hit_evalue.get((qp, pid, sp), '')
            tgt_row[sp] = '' if sp == qp else hit_target.get((qp, pid, sp), '')
        rows.append(row)
        evalue_rows.append(ev_row)
        target_rows.append(tgt_row)

    matrix = pd.DataFrame(rows, columns=columns)
    matrix.to_csv(args.output_matrix, sep='\t', index=False)

    if args.output_evalues:
        evalues_df = pd.DataFrame(evalue_rows, columns=columns)
        evalues_df.to_csv(args.output_evalues, sep='\t', index=False)

    if args.output_targets:
        targets_df = pd.DataFrame(target_rows, columns=columns)
        targets_df.to_csv(args.output_targets, sep='\t', index=False)

    n_query     = len(query_ids)
    query_cols  = sorted(query_ids)
    other_cols  = sorted(other_ids)
    if matrix.empty:
        candidates = []
    else:
        query_count = matrix[query_cols].sum(axis=1)
        other_count = matrix[other_cols].sum(axis=1) if other_cols else 0
        other_frac = (other_count / len(other_cols)) if other_cols else 0
        keep = (query_count / n_query >= args.ingroup_min_frac) & (other_frac <= args.other_max_frac)
        kept = matrix[keep]
        candidates = (kept['source_proteome'] + '::' + kept['protein_id']).tolist()

    with open(args.output_candidates, 'w') as fh:
        if candidates:
            fh.write('\n'.join(candidates) + '\n')

    direction = 'ingroup' if args.query_group == 'IN' else 'outgroup'
    print(f"Candidates: {len(candidates)} / {len(matrix)} {direction} proteins pass thresholds",
          file=sys.stderr)


if __name__ == '__main__':
    main()
