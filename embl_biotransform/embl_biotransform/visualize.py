"""
visualize.py
============

Two complementary views of the same evidence:

- `plot_evidence_graph`  -- the actual network (sample -> taxa/functions ->
  reference statements -> reaction classes), colour-coded by node type, so
  you can see *how* a prediction was reached, not just its score.
- `plot_evidence_scores` -- a horizontal bar chart of predicted reaction
  classes ranked by confidence, annotated with how many Tier-1 (enzyme)
  vs Tier-2 (microbe-only) evidence items support each one.

Both take the `EvidenceGraph` / predictions directly and save a PNG (and
return the Figure, so they also work inline in a notebook).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx

NODE_COLORS = {
    "sample": "#4C72B0",
    "taxon": "#55A868",
    "functional_annotation": "#C44E52",
    "reference": "#8172B2",
    "reaction_class": "#CCB974",
    "compound": "#64B5CD",
}


def plot_evidence_graph(evidence_graph, sample_node: str | None = None,
                         save_path: str | None = None, figsize=(11, 8)):
    """Draws the evidence graph. If `sample_node` is given, restricts the
    view to that sample's connected component so a multi-sample graph
    doesn't turn into an unreadable hairball."""
    g = evidence_graph.graph
    if sample_node is not None:
        component = nx.node_connected_component(g.to_undirected(), sample_node)
        g = g.subgraph(component)

    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(g, seed=7, k=0.9)

    for ntype, color in NODE_COLORS.items():
        nodes = [n for n, d in g.nodes(data=True) if d.get("type") == ntype]
        if not nodes:
            continue
        nx.draw_networkx_nodes(g, pos, nodelist=nodes, node_color=color,
                                node_size=550, label=ntype, ax=ax, alpha=0.9)

    # colour Tier-1 vs Tier-2 reference nodes differently (outline) so the
    # evidence strength is visible directly on the graph
    tier1 = [n for n, d in g.nodes(data=True) if d.get("type") == "reference" and d.get("tier") == 1]
    tier2 = [n for n, d in g.nodes(data=True) if d.get("type") == "reference" and d.get("tier") == 2]
    if tier1:
        nx.draw_networkx_nodes(g, pos, nodelist=tier1, node_color="none",
                                edgecolors="#2E2530", linewidths=2.2, node_size=650, ax=ax)
    if tier2:
        coll = nx.draw_networkx_nodes(g, pos, nodelist=tier2, node_color="none",
                                       edgecolors="#2E2530", linewidths=1.4, node_size=650, ax=ax)
        coll.set_linestyle((0, (3, 2)))  # dashed outline = Tier 2 (microbe-only)

    nx.draw_networkx_edges(g, pos, alpha=0.35, arrows=True, arrowsize=10, ax=ax)
    labels = {n: (d.get("organism") or d.get("description") or d.get("entity_name")
                  or d.get("name") or d.get("sample_id") or n.split(":")[-1])
              for n, d in g.nodes(data=True)}
    nx.draw_networkx_labels(g, pos, labels=labels, font_size=7, ax=ax)

    ax.set_title("Biotransformation evidence graph\n"
                 "(solid outline = Tier 1 enzyme evidence, dashed = Tier 2 microbe-only evidence)")
    ax.legend(scatterpoints=1, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    ax.axis("off")
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_evidence_scores(predictions, save_path: str | None = None, figsize=(9, 5)):
    """Ranked bar chart of predicted reaction classes with Tier-1/Tier-2
    evidence counts annotated on each bar."""
    predictions = list(predictions)[::-1]  # so the top prediction plots at the top
    labels = [p.reaction_class for p in predictions]
    scores = [p.score for p in predictions]
    colors = ["#8172B2" if p.tier1_hits else "#CCB974" for p in predictions]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.barh(labels, scores, color=colors)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Confidence score")
    ax.set_title("Predicted biotransformations, ranked by evidence-backed confidence")

    for bar, p in zip(bars, predictions):
        note = f"{len(p.tier1_hits)} enzyme hit(s), {len(p.tier2_hits)} microbe hit(s)"
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                 note, va="center", fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig
