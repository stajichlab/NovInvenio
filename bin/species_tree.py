#!/usr/bin/env python3
"""Root a species tree (Newick) on the OUT group and draw it as inline SVG.

Rooting: OUT means GROUP == OUT only (other roles such as NEAR_INGROUP stay
ordinary tips). Pick the edge whose split best separates OUT tips from all other tips
(score = OUT tips on one side minus non-OUT tips on that side; ties go to the
longest edge) and root at the middle of that edge. If the tree has no OUT tips,
midpoint-root. `n_out_outside` counts OUT tips that are not on the OUT side of
the chosen split.

Tip labels (samplesheet `Short`) may contain `.`, `:`, `|` and spaces. The
parser reads quoted labels and writes quotes when a label needs them.

Usage:
  species_tree.py --newick tree.nwk --config config.csv --out_prefix species_tree
Option --tip_labels species draws the samplesheet Species name instead of Short
(duplicates get ` <Strain>` or ` (<Short>)`). Only the SVG text changes.
Writes <prefix>.rooted.nwk, <prefix>.svg and <prefix>.rooting.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from html import escape
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "lib"))

_SPECIAL = set(" \t\n,:;()[]'")


class Node:
    def __init__(self, label: str = "", length: float = 0.0):
        self.label = label
        self.length = length
        self.children: list[Node] = []


def parse_newick(text: str) -> Node:
    """Parse one Newick tree. Handles quoted labels and missing lengths (0)."""
    s = text.strip()
    pos = 0
    n = len(s)

    def read_label() -> str:
        nonlocal pos
        if pos < n and s[pos] == "'":
            pos += 1
            out = []
            while pos < n:
                if s[pos] == "'":
                    if pos + 1 < n and s[pos + 1] == "'":
                        out.append("'")
                        pos += 2
                        continue
                    pos += 1
                    break
                out.append(s[pos])
                pos += 1
            return "".join(out)
        start = pos
        while pos < n and s[pos] not in ",:;()":
            pos += 1
        return s[start:pos].strip()

    def read_length() -> float:
        nonlocal pos
        if pos < n and s[pos] == ":":
            pos += 1
            start = pos
            while pos < n and s[pos] not in ",;()":
                pos += 1
            tok = s[start:pos].strip()
            return float(tok) if tok else 0.0
        return 0.0

    root = Node()
    stack: list[Node] = []
    cur = root
    while pos < n:
        c = s[pos]
        if c == "(":
            child = Node()
            cur.children.append(child)
            stack.append(cur)
            cur = child
            pos += 1
        elif c == ",":
            parent = stack[-1]
            child = Node()
            parent.children.append(child)
            cur = child
            pos += 1
        elif c == ")":
            cur = stack.pop()
            pos += 1
            cur.label = read_label()
            cur.length = read_length()
        elif c == ";":
            break
        elif c.isspace():
            pos += 1
        else:
            cur.label = read_label()
            cur.length = read_length()
    return root


def _fmt_label(label: str) -> str:
    if any(ch in _SPECIAL for ch in label):
        return "'" + label.replace("'", "''") + "'"
    return label


def write_newick(node: Node, top: bool = True) -> str:
    if node.children:
        body = "(" + ",".join(write_newick(c, False) for c in node.children) + ")"
    else:
        body = ""
    body += _fmt_label(node.label) if node.label else ""
    if not top:
        body += f":{node.length:.10g}"
    return body + (";" if top else "")


def tip_labels(node: Node) -> list[str]:
    out, stack = [], [node]
    while stack:
        x = stack.pop()
        if x.children:
            stack.extend(reversed(x.children))
        else:
            out.append(x.label)
    return out


def tip_distance(root: Node, a: str, b: str) -> float:
    """Path length between two tips (ignores the root's own length)."""
    paths: dict[str, list[tuple[Node, float]]] = {}

    def walk(node: Node, trail: list[Node]):
        trail = trail + [node]
        if not node.children:
            paths[node.label] = trail
        for c in node.children:
            walk(c, trail)

    walk(root, [])
    pa, pb = paths[a], paths[b]
    i = 0
    while i < min(len(pa), len(pb)) and pa[i] is pb[i]:
        i += 1
    return sum(x.length for x in pa[i:]) + sum(x.length for x in pb[i:])


# ---- unrooted graph view ------------------------------------------------

def _to_graph(root: Node):
    """Return (adj, labels, edges). Nodes are ints; a binary root is suppressed."""
    adj: dict[int, list[tuple[int, float]]] = {}
    labels: dict[int, str] = {}
    counter = [0]

    def new(label: str) -> int:
        i = counter[0]
        counter[0] += 1
        adj[i] = []
        labels[i] = label
        return i

    def link(a: int, b: int, w: float):
        adj[a].append((b, w))
        adj[b].append((a, w))

    def build(node: Node, parent: int | None, w: float):
        i = new(node.label if not node.children else "")
        if parent is not None:
            link(parent, i, w)
        for c in node.children:
            build(c, i, c.length)
        return i

    if len(root.children) == 2:
        left, right = root.children
        a = build(left, None, 0.0)
        b = build(right, None, 0.0)
        link(a, b, left.length + right.length)
    else:
        build(root, None, 0.0)
    edges = []
    for u, nbrs in adj.items():
        for v, w in nbrs:
            if u < v:
                edges.append((u, v, w))
    return adj, labels, edges


def _side(adj, start: int, block: int) -> set[int]:
    seen, stack = {start}, [start]
    while stack:
        x = stack.pop()
        for y, _ in adj[x]:
            if y != block and y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


def _reroot(adj, labels, u: int, v: int, w: float, dist_u: float | None = None) -> Node:
    """Root on edge u-v at distance dist_u from u (default: middle)."""
    du = w / 2 if dist_u is None else dist_u

    def build(x: int, parent: int, length: float) -> Node:
        node = Node(labels[x], length)
        for y, ww in adj[x]:
            if y != parent:
                node.children.append(build(y, x, ww))
        return node

    root = Node()
    root.children = [build(u, v, du), build(v, u, w - du)]
    return root


def _is_tip(adj, i: int) -> bool:
    return len(adj[i]) == 1


def _midpoint_edge(adj, labels, edges):
    tips = [i for i in adj if _is_tip(adj, i)]

    def far(src: int):
        dist = {src: 0.0}
        prev = {}
        stack = [src]
        while stack:
            x = stack.pop()
            for y, ww in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + ww
                    prev[y] = x
                    stack.append(y)
        best = max(tips, key=lambda t: dist[t])
        return best, dist, prev

    a, _, _ = far(tips[0])
    b, dist, prev = far(a)
    half = dist[b] / 2
    # walk from b towards a until the cumulative distance from a passes half
    x = b
    while x != a:
        p = prev[x]
        if dist[p] <= half <= dist[x]:
            return p, x, dist[x] - dist[p], half - dist[p]
        x = p
    return None


def root_on_outgroup(newick: str, groups: dict[str, str]) -> tuple[str, dict]:
    """Root `newick`; `groups` maps tip label -> GROUP. Returns (newick, info)."""
    tree = parse_newick(newick)
    adj, labels, edges = _to_graph(tree)
    tip_ids = [i for i in adj if _is_tip(adj, i)]
    is_out = {i: groups.get(labels[i]) == "OUT" for i in tip_ids}
    n_out = sum(is_out.values())
    n_tips = len(tip_ids)

    if n_tips < 2 or not edges:
        info = {"method": "midpoint", "n_out": n_out, "n_out_outside": n_out, "n_non_out_inside": 0}
        return write_newick(tree), info

    if n_out == 0:
        res = _midpoint_edge(adj, labels, edges)
        u, v, w, du = res
        rooted = _reroot(adj, labels, u, v, w, du)
        return write_newick(rooted), {"method": "midpoint", "n_out": 0, "n_out_outside": 0, "n_non_out_inside": 0}

    best = None
    for u, v, w in edges:
        side_v = _side(adj, v, u)
        tips_v = [t for t in side_v if t in is_out]
        out_v = sum(is_out[t] for t in tips_v)
        non_v = len(tips_v) - out_v
        score_v = out_v - non_v
        score_u = (n_out - out_v) - ((n_tips - n_out) - non_v)
        # orient so the OUT side comes first (the u side of the new root)
        if score_v >= score_u:
            cand = (score_v, w, v, u, out_v, non_v)
        else:
            cand = (score_u, w, u, v, n_out - out_v, (n_tips - n_out) - non_v)
        if best is None or cand[:2] > best[:2]:
            best = cand
    _, w, out_node, other, out_in, non_in = best
    rooted = _reroot(adj, labels, out_node, other, w)
    pure = out_in == n_out and non_in == 0
    info = {
        "method": "outgroup" if pure else "outgroup_best_split",
        "n_out": n_out,
        "n_out_outside": n_out - out_in,
        "n_non_out_inside": non_in,
    }
    return write_newick(rooted), info


# ---- SVG ------------------------------------------------------------------

_PALETTE = {
    "IN": "#1f77b4",
    "OUT": "#d95f02",
    "DISCOVERY_TARGET": "#1f77b4",
    "DISCOVERY_OUT": "#d95f02",
    "NEAR_INGROUP": "#7570b3",
    "BROAD_OUTGROUP": "#a6761d",
}
_FALLBACK = ["#1b9e77", "#e7298a", "#66a61e", "#666666"]
_ROW = 16
_FONT = 12
_CHAR_W = 7.0


def render_svg(newick: str, groups: dict[str, str],
               labels: dict[str, str] | None = None) -> str:
    """Draw the tree. `labels` maps tip ID -> display text; a missing tip shows its ID."""
    root = parse_newick(newick)
    tips = tip_labels(root)
    labels = labels or {}

    def shown(tip: str) -> str:
        return labels.get(tip) or tip

    colour: dict[str, str] = {}
    for g in sorted(set(groups.get(t, "unknown") for t in tips)):
        colour[g] = _PALETTE.get(g) or _FALLBACK[len(
            [k for k in colour if k not in _PALETTE]) % len(_FALLBACK)]

    xs: dict[int, float] = {}
    ys: dict[int, float] = {}
    order: list[Node] = []
    counter = [0]

    def lay(node: Node, depth: float):
        xs[id(node)] = depth
        order.append(node)
        if node.children:
            for c in node.children:
                lay(c, depth + c.length)
            ys[id(node)] = (ys[id(node.children[0])] + ys[id(node.children[-1])]) / 2
        else:
            ys[id(node)] = counter[0] * _ROW + _ROW
            counter[0] += 1

    lay(root, 0.0)
    max_depth = max(xs.values()) or 1.0
    plot_w = 500.0
    left = 20.0
    label_w = max(len(shown(t)) for t in tips) * _CHAR_W + 20
    legend_h = 24 + 18 * len(colour)
    height = counter[0] * _ROW + _ROW + 20
    height = max(height, legend_h + 10)
    width = left + plot_w + label_w + 160
    sx = plot_w / max_depth

    def px(node: Node) -> float:
        return left + xs[id(node)] * sx

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" '
        f'viewBox="0 0 {width:.0f} {height:.0f}" font-family="sans-serif" font-size="{_FONT}">',
        f'<rect width="100%" height="100%" fill="#ffffff"/>',
        '<g stroke="#444444" stroke-width="1" fill="none">',
    ]
    for node in order:
        if node.children:
            y1, y2 = ys[id(node.children[0])], ys[id(node.children[-1])]
            parts.append(f'<line x1="{px(node):.1f}" y1="{y1:.1f}" x2="{px(node):.1f}" y2="{y2:.1f}"/>')
            for c in node.children:
                parts.append(
                    f'<line x1="{px(node):.1f}" y1="{ys[id(c)]:.1f}" x2="{px(c):.1f}" y2="{ys[id(c)]:.1f}"/>')
    parts.append("</g>")
    for node in order:
        if not node.children:
            g = groups.get(node.label, "unknown")
            parts.append(
                f'<text class="tip" x="{px(node) + 6:.1f}" y="{ys[id(node)] + 4:.1f}" '
                f'fill="{colour[g]}">{escape(shown(node.label), quote=False)}</text>')
    lx = left + plot_w + label_w + 10
    parts.append(f'<g class="legend"><text x="{lx:.0f}" y="16" font-weight="bold">GROUP</text>')
    for i, (g, c) in enumerate(colour.items()):
        y = 34 + 18 * i
        parts.append(f'<rect x="{lx:.0f}" y="{y - 10}" width="12" height="12" fill="{c}"/>')
        parts.append(f'<text x="{lx + 18:.0f}" y="{y}">{escape(g, quote=False)}</text>')
    parts.append("</g></svg>")
    return "\n".join(parts) + "\n"


def species_labels(samples) -> dict[str, str]:
    """Short -> unique display name from the Species column.

    Tips that share a Species value get ` <Strain>` (or ` (<Short>)` when Strain is
    empty). Empty Species gives no entry, so the tip keeps its Short ID.
    """
    counts: dict[str, int] = {}
    for s in samples:
        if s.species:
            counts[s.species] = counts.get(s.species, 0) + 1
    out: dict[str, str] = {}
    for s in samples:
        if not s.species:
            continue
        if counts[s.species] > 1:
            out[s.short] = f"{s.species} {s.strain}" if s.strain else f"{s.species} ({s.short})"
        else:
            out[s.short] = s.species
    return out


def main(argv=None) -> int:
    from config_parser import parse_config

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--newick", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out_prefix", default="species_tree")
    ap.add_argument("--tip_labels", choices=("short", "species"), default="short",
                    help="SVG tip text: samplesheet Short ID (default) or Species name")
    a = ap.parse_args(argv)

    samples = parse_config(a.config)
    groups = {s.short: s.group for s in samples}
    labels = species_labels(samples) if a.tip_labels == "species" else None
    text = Path(a.newick).read_text()
    if "OUT" in groups.values() and not set(tip_labels(parse_newick(text))) & set(groups):
        print("tree tip labels do not match config Short names", file=sys.stderr)
        return 3
    rooted, info = root_on_outgroup(text, groups)
    Path(f"{a.out_prefix}.rooted.nwk").write_text(rooted + "\n")
    Path(f"{a.out_prefix}.svg").write_text(render_svg(rooted, groups, labels))
    Path(f"{a.out_prefix}.rooting.json").write_text(json.dumps(info, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
