
'''
The cornish plotting stack (SPEC-09A).

``SkyPlot`` drives ``Ast.Plot`` through a recording Grf backend and captures the
result as a replayable draw-stream document; SVG, JSON, the browser and
matplotlib are all consumers of that one stream (D17).

Nothing here imports matplotlib or any transport library at module level —
``render_svg`` is a pure call, which is what makes the recommended production
shape (Flask pattern 1) a one-liner. Both properties are gate-enforced by
``test_render_svg.py``.

``cornish.plot.matplotlib`` still holds the frozen v1 ``SkyPlot``, superseded by
this one and deliberately NOT re-exported here: it imports matplotlib at module
level, and importing it would break the "renders with no matplotlib installed"
guarantee. It is removed at 2.0.
'''

from .envelope import (SCHEMA_VERSION, GrfOp, LayerAction, PrimType, canonical_json,
                       envelope_hash, make_envelope, validate_document, validate_envelope)
from .metrics import FontMetrics
from .palette import DEFAULT_PALETTE, normalize_colour, normalize_palette_entry
from .serial import GrfSequence
from .sky import SkyPlot
from .svg import SVGRenderer, render_svg

__all__ = [
	'SkyPlot', 'GrfSequence', 'render_svg', 'SVGRenderer', 'FontMetrics',
	'DEFAULT_PALETTE', 'normalize_colour', 'normalize_palette_entry',
	'SCHEMA_VERSION', 'GrfOp', 'LayerAction', 'PrimType',
	'canonical_json', 'envelope_hash', 'make_envelope',
	'validate_document', 'validate_envelope',
]
