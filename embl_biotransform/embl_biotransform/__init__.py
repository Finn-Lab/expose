"""
embl_biotransform
==================

A small, extensible toolkit that pulls compound / reaction / enzyme
evidence from EMBL-EBI resources (ChEBI, Rhea, UniProt, ChEMBL),
integrates it into a single evidence graph, scores candidate
biotransformations for a query compound, and renders the evidence
behind each prediction.

See README.md for the architecture and how to extend it.
"""

from .fetchers import ChEBIClient, RheaClient, UniProtClient, ChEMBLClient, MGnifyClient
from .mgnipy_client import MGnipySource
from .reference import load_reference_associations, parse_reference_text, ReferenceFormatError
from .integrate import EvidenceGraph
from .predict import BiotransformationPredictor, BiotransformationPrediction
from .pipeline import (
    run_analysis, run_multi_analysis, build_grid,
    parse_detections_text, DetectionFormatError,
)
from .biosiftr import (
    BioSIFTRFormatError, load_run as load_biosiftr_run,
    parse_table as parse_biosiftr_table, summarise_run as summarise_biosiftr_run,
)
from .chembl_references import (
    ChEMBLBiotransformationSource, ChEMBLExtraction, ChEMBLReferenceError,
    build_reference_associations, classify_assay,
)
from .visualize import plot_evidence_graph, plot_evidence_scores

__all__ = [
    "ChEBIClient",
    "RheaClient",
    "UniProtClient",
    "ChEMBLClient",
    "MGnifyClient",
    "MGnipySource",
    "load_reference_associations",
    "parse_reference_text",
    "ReferenceFormatError",
    "EvidenceGraph",
    "BiotransformationPredictor",
    "BiotransformationPrediction",
    "run_analysis",
    "run_multi_analysis",
    "build_grid",
    "parse_detections_text",
    "DetectionFormatError",
    "BioSIFTRFormatError",
    "load_biosiftr_run",
    "parse_biosiftr_table",
    "summarise_biosiftr_run",
    "ChEMBLBiotransformationSource",
    "ChEMBLExtraction",
    "ChEMBLReferenceError",
    "build_reference_associations",
    "classify_assay",
    "plot_evidence_graph",
    "plot_evidence_scores",
]

__version__ = "0.3.0"
