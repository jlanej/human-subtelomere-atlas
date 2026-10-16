"""FracMinHash sketches and Mash-style ANI. Implementation: atlas/subtel_atlas/sketch.py (shared with the query tool)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'atlas'))
from subtel_atlas.sketch import sketch, ani_from_jaccard, _splitmix, _LUT  # noqa: E402,F401
