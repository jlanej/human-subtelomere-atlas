"""Telomere block caller (tip frame; (CCCTAA)n). Implementation: atlas/subtel_atlas/telomere.py (shared with the query tool)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'atlas'))
from subtel_atlas.telomere import CAN, find_block, _rots, _covmask, _frac  # noqa: E402,F401
