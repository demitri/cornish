
'''
The Grf checklist as executable spec (SPEC-09A §7.1).

AST talks to a graphics backend through a fixed set of callbacks. This module
runs that contract rather than describing it, and — where the shared subset
allows — runs it against ``starlink.Grf.grf_matplotlib`` too, so the tests
document what AST actually requires rather than what cornish happens to
implement.
'''

import math

import pytest
import starlink.Ast as Ast

from conftest import golden_ast_plot
from cornish.plot.envelope import GrfOp
from cornish.plot.palette import DEFAULT_PALETTE
from cornish.plot.serial import GrfSequence


#: The full Grf contract AST requires, plus pyast's two colour hooks.
GRF_METHODS = ("Line", "Mark", "Text", "TxExt", "Qch", "Scales", "Attr", "Cap",
               "BBuf", "EBuf", "Flush", "ColToInt", "IntToCol")

ATTRIBUTE_CODES = (Ast.grfSTYLE, Ast.grfWIDTH, Ast.grfSIZE, Ast.grfFONT, Ast.grfCOLOUR)
PRIMITIVE_CODES = (Ast.grfLINE, Ast.grfMARK, Ast.grfTEXT)


#: One shared reference backend. Every use below is read-only (surface checks,
#: Cap/Qch/TxExt queries), and a pyplot figure per test would leak figures past
#: matplotlib's open-figure warning threshold — which the v1 SkyPlot tests are
#: already near.
_MATPLOTLIB_BACKEND = []


def matplotlib_backend():
	'''
	A ``grf_matplotlib`` on a headless Agg figure; skips when matplotlib is
	unavailable. Used only for the shared-subset tests.
	'''
	matplotlib = pytest.importorskip("matplotlib")
	if not _MATPLOTLIB_BACKEND:
		matplotlib.use("Agg", force=True)
		# pyplot (not a bare Figure) so the figure carries a FigureCanvasAgg:
		# grf_matplotlib measures text by asking the canvas for a renderer,
		# which is precisely the renderer-dependence GrfSequence avoids (D17).
		from matplotlib import pyplot
		from starlink.Grf import grf_matplotlib
		axes = pyplot.figure().add_axes([0, 0, 1, 1])
		_MATPLOTLIB_BACKEND.append(grf_matplotlib(axes))
	return _MATPLOTLIB_BACKEND[0]


@pytest.fixture(params=["GrfSequence", "grf_matplotlib"])
def any_backend(request):
	''' Both backends, for the parts of the contract they genuinely share. '''
	if request.param == "GrfSequence":
		return GrfSequence()
	return matplotlib_backend()


# ---------------------------------------------------------------------------
# the contract surface
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", GRF_METHODS)
def test_every_grf_method_exists(name):
	assert callable(getattr(GrfSequence(), name)), f"GrfSequence is missing the Grf method {name}()"


@pytest.mark.parametrize("name", GRF_METHODS)
def test_grf_matplotlib_shares_the_same_method_surface(name):
	''' The contract is AST's, not cornish's — the reference backend has it too. '''
	assert callable(getattr(matplotlib_backend(), name))


def test_buffering_hooks_are_callable_and_return_none():
	backend = GrfSequence()
	assert backend.BBuf() is None
	assert backend.EBuf() is None
	assert backend.Flush() is None


# ---------------------------------------------------------------------------
# Attr: the query protocol
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("attribute", ATTRIBUTE_CODES)
@pytest.mark.parametrize("primitive", PRIMITIVE_CODES)
def test_attr_query_returns_current_value_without_mutating(attribute, primitive):
	'''
	``Attr(attr, Ast.BAD, prim)`` is AST's *query* form: report, do not change.
	'''
	backend = GrfSequence()
	backend.Attr(attribute, 3, primitive)
	first = backend.Attr(attribute, Ast.BAD, primitive)
	second = backend.Attr(attribute, Ast.BAD, primitive)
	assert first == second == 3.0


@pytest.mark.parametrize("attribute", ATTRIBUTE_CODES)
@pytest.mark.parametrize("primitive", PRIMITIVE_CODES)
def test_attr_set_query_restore_round_trips(attribute, primitive):
	''' Every (attr, prim) pair round-trips — AST relies on save/restore. '''
	backend = GrfSequence()
	original = backend.Attr(attribute, Ast.BAD, primitive)
	previous = backend.Attr(attribute, 4, primitive)
	assert previous == original
	assert backend.Attr(attribute, Ast.BAD, primitive) == 4.0
	backend.Attr(attribute, original, primitive)
	assert backend.Attr(attribute, Ast.BAD, primitive) == original


def test_attr_setting_one_primitive_does_not_disturb_another():
	backend = GrfSequence()
	backend.Attr(Ast.grfCOLOUR, 9, Ast.grfLINE)
	assert backend.Attr(Ast.grfCOLOUR, Ast.BAD, Ast.grfTEXT) == 1.0


def test_attr_rejects_unknown_codes():
	''' An unknown code is a wiring bug; guessing one would hide it. '''
	backend = GrfSequence()
	with pytest.raises(ValueError):
		backend.Attr(99, 1, Ast.grfLINE)
	with pytest.raises(ValueError):
		backend.Attr(Ast.grfCOLOUR, 1, 99)


# ---------------------------------------------------------------------------
# Cap / Qch / Scales
# ---------------------------------------------------------------------------

def test_cap_values(any_backend):
	'''
	``ESC = 0`` is load-bearing: it is what makes AST deliver label text plain
	instead of sprinkled with ``%`` escape sequences (§3.3).
	'''
	assert any_backend.Cap(Ast.grfSCALES, 0) == 1
	assert any_backend.Cap(Ast.grfMJUST, 0) == 1
	assert any_backend.Cap(Ast.grfESC, 1) == 0


def test_qch_returns_two_positive_floats(any_backend):
	height = any_backend.Qch()
	assert len(height) == 2
	assert all(isinstance(value, float) and value > 0 for value in height)


def test_qch_is_isotropic_for_the_recorder():
	''' The gbox is isotropic by construction (§2.2), so chv == chh. '''
	backend = GrfSequence(char_height=0.03)
	assert backend.Qch() == (0.03, 0.03)


def test_scales_is_unity():
	''' (1, 1) is an honest answer only because the gbox is isotropic (§2.2). '''
	assert GrfSequence().Scales() == (1.0, 1.0)


def test_char_height_must_be_positive():
	for bad in (0, -1, float("nan")):
		with pytest.raises(ValueError):
			GrfSequence(char_height=bad)


# ---------------------------------------------------------------------------
# TxExt
# ---------------------------------------------------------------------------

def test_txext_returns_eight_finite_floats(any_backend):
	extents = any_backend.TxExt("label", 0.5, 0.5, "CC", 0.0, 1.0)
	assert len(extents) == 8
	assert all(math.isfinite(value) for value in extents)


def test_txext_width_grows_monotonically_with_string_length():
	backend = GrfSequence()
	widths = []
	for text in ("i", "ii", "iiii", "iiiiiiii"):
		extents = backend.TxExt(text, 0.0, 0.0, "BL", 0.0, 1.0)
		widths.append(max(extents[:4]) - min(extents[:4]))
	assert widths == sorted(widths)
	assert widths[0] < widths[-1]


def test_txext_of_an_empty_string_is_degenerate_at_the_anchor():
	extents = GrfSequence().TxExt("", 0.25, 0.75, "CC", 0.0, 1.0)
	assert extents == (0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75)


def test_txext_scales_with_the_text_size_attribute():
	''' The size attr multiplies the em, so TxExt must consult it. '''
	backend = GrfSequence()
	small = backend.TxExt("AAA", 0.0, 0.0, "BL", 0.0, 1.0)
	backend.Attr(Ast.grfSIZE, 2.0, Ast.grfTEXT)
	large = backend.TxExt("AAA", 0.0, 0.0, "BL", 0.0, 1.0)
	small_width = max(small[:4]) - min(small[:4])
	large_width = max(large[:4]) - min(large[:4])
	assert large_width == pytest.approx(2 * small_width)


def test_txext_rejects_an_illegal_justification():
	with pytest.raises(ValueError):
		GrfSequence().TxExt("x", 0.0, 0.0, "ZZ", 0.0, 1.0)


# ---------------------------------------------------------------------------
# ColToInt / IntToCol
# ---------------------------------------------------------------------------

def test_coltoint_passes_integers_through():
	assert GrfSequence().ColToInt(3) == 3


def test_coltoint_grows_the_palette_and_round_trips():
	backend = GrfSequence()
	index = backend.ColToInt("#4a7f7b")
	assert index == len(DEFAULT_PALETTE)          # appended past the base palette
	assert backend.IntToCol(index) == "#4a7f7b"
	assert backend.ColToInt("#4a7f7b") == index   # idempotent, not appended twice


def test_coltoint_resolves_css_named_colours():
	backend = GrfSequence()
	assert backend.IntToCol(backend.ColToInt("red")) == "#ff0000"


def test_coltoint_returns_minus_one_for_garbage():
	'''
	−1 is the Grf contract's "not a legal colour" — AST reports the error
	itself, so nothing is silently substituted.
	'''
	assert GrfSequence().ColToInt("not-a-colour") == -1


def test_inttocol_out_of_range_is_none():
	assert GrfSequence().IntToCol(999) is None


def test_palette_property_is_a_copy():
	backend = GrfSequence()
	snapshot = backend.palette
	snapshot.append("#123456")
	assert backend.palette == list(DEFAULT_PALETTE)


def test_colour_attribute_is_recorded_as_an_integer_index():
	'''
	AST hands every Grf attribute over as a double, but colour is a palette
	INDEX — recording 11.0 would violate the schema and defeat the
	index-vs-CSS-string resolution done at take time.
	'''
	backend = GrfSequence()
	backend.Attr(Ast.grfCOLOUR, 11.0, Ast.grfLINE)
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	attribute_command = backend.take_commands()[0]
	assert attribute_command["colour"] == 11
	assert isinstance(attribute_command["colour"], int)


# ---------------------------------------------------------------------------
# attribute coalescing and take_commands
# ---------------------------------------------------------------------------

def test_attr_chatter_is_coalesced_into_one_net_command():
	'''
	AST sets and restores attributes constantly. Only the net state that
	actually differs, immediately before a draw op, reaches the stream.
	'''
	backend = GrfSequence()
	for value in (2, 3, 4, 5):                 # four sets, no draw in between
		backend.Attr(Ast.grfCOLOUR, value, Ast.grfLINE)
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	commands = backend.take_commands()
	assert [command["op"] for command in commands] == [GrfOp.ATTR.value, GrfOp.LINE.value]
	assert commands[0]["colour"] == 5          # the net value, not the four steps


def test_no_attr_command_when_nothing_changed():
	backend = GrfSequence()
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	backend.Line(2, [0.0, 1.0], [1.0, 0.0])
	assert [command["op"] for command in backend.take_commands()] == [GrfOp.LINE.value] * 2


def test_attr_set_then_restored_before_any_draw_emits_nothing():
	''' A save/restore pair that brackets no drawing is pure chatter. '''
	backend = GrfSequence()
	original = backend.Attr(Ast.grfCOLOUR, 7, Ast.grfLINE)
	backend.Attr(Ast.grfCOLOUR, original, Ast.grfLINE)
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	assert [command["op"] for command in backend.take_commands()] == [GrfOp.LINE.value]


def test_take_commands_resets_the_attr_baseline():
	'''
	§3.4/D18: each taken chunk becomes one layer and must be self-contained, so
	non-default state is re-emitted at the head of the next chunk. Otherwise
	deleting one layer would restyle the next.
	'''
	backend = GrfSequence()
	backend.Attr(Ast.grfCOLOUR, 6, Ast.grfLINE)
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	backend.take_commands()

	backend.Line(2, [0.0, 1.0], [1.0, 0.0])    # same state, new chunk
	second = backend.take_commands()
	assert second[0]["op"] == GrfOp.ATTR.value
	assert second[0]["colour"] == 6


def test_take_commands_pops():
	backend = GrfSequence()
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	assert len(backend.take_commands()) == 1
	assert backend.take_commands() == []


def test_grown_palette_colours_serialize_as_literal_css():
	'''
	D18: a base-palette index stays an int so consumers can re-theme it; an
	index grown by ColToInt is written out as its CSS string, so the layer does
	not depend on a palette the consumer never received.
	'''
	backend = GrfSequence()
	grown = backend.ColToInt("#4a7f7b")
	backend.Attr(Ast.grfCOLOUR, grown, Ast.grfLINE)
	backend.Line(2, [0.0, 1.0], [0.0, 1.0])
	assert backend.take_commands()[0]["colour"] == "#4a7f7b"


def test_base_palette_colours_stay_integer_indices():
	backend = GrfSequence()
	backend.Attr(Ast.grfCOLOUR, 4, Ast.grfMARK)
	backend.Mark(1, [0.5], [0.5], 4)
	assert backend.take_commands()[0]["colour"] == 4


# ---------------------------------------------------------------------------
# recording the primitives
# ---------------------------------------------------------------------------

def test_line_records_coordinates_verbatim():
	backend = GrfSequence()
	backend.Line(3, [0.0, 0.5, 1.0], [1.0, 0.5, 0.0])
	command = backend.take_commands()[0]
	assert command == {"op": "line", "x": [0.0, 0.5, 1.0], "y": [1.0, 0.5, 0.0]}


def test_line_honours_the_point_count():
	''' AST passes buffers that may be longer than n. '''
	backend = GrfSequence()
	backend.Line(2, [0.0, 1.0, 9.0], [0.0, 1.0, 9.0])
	assert backend.take_commands()[0]["x"] == [0.0, 1.0]


def test_degenerate_primitives_draw_nothing():
	backend = GrfSequence()
	backend.Line(1, [0.0], [0.0])
	backend.Mark(0, [], [], 4)
	backend.Text("", 0.0, 0.0, "CC", 0.0, 1.0)
	assert backend.take_commands() == []


def test_bad_coordinates_are_refused_loudly():
	'''
	``Ast.BAD`` is a *finite* float (-1.797e308), so it would pass any ordinary
	finiteness check and be recorded as a real coordinate. AST breaks its own
	polylines at bad points and should never hand one over; if it ever does,
	that must surface rather than deform the geometry silently.
	'''
	backend = GrfSequence()
	with pytest.raises(ValueError, match="Ast.BAD"):
		backend.Line(2, [0.0, Ast.BAD], [0.0, 1.0])
	with pytest.raises(ValueError):
		backend.Mark(1, [float("nan")], [0.0], 4)


def test_text_records_the_full_anchor_specification():
	backend = GrfSequence()
	backend.Text("M51", 0.25, 0.75, "bc", 0.0, 1.0)
	command = backend.take_commands()[0]
	assert command == {"op": "text", "text": "M51", "x": 0.25, "y": 0.75,
	                   "just": "BC", "up": [0.0, 1.0]}


def test_text_rejects_an_illegal_justification():
	with pytest.raises(ValueError):
		GrfSequence().Text("x", 0.0, 0.0, "QQ", 0.0, 1.0)


# ---------------------------------------------------------------------------
# the smoke test: pyast actually accepts the object
# ---------------------------------------------------------------------------

def test_ast_plot_grid_drives_the_recorder():
	'''
	The test that matters most: a real ``Ast.Plot`` drawing a real grid into
	``GrfSequence`` must raise nothing and produce both geometry and labels.
	'''
	backend = GrfSequence()
	plot = golden_ast_plot(backend)
	plot.grid()
	commands = backend.take_commands()

	lines = [command for command in commands if command["op"] == GrfOp.LINE.value]
	texts = [command for command in commands if command["op"] == GrfOp.TEXT.value]
	assert lines, "a grid produced no line commands"
	assert texts, "a grid produced no text commands"


def test_grid_labels_arrive_without_escape_sequences():
	'''
	The consequence of ``Cap(grfESC) = 0`` [verified in the design session:
	0 of 17 labels contained a ``%`` escape]. Envelope consumers MUST NOT
	interpret ``%`` sequences, so producers must not emit them.
	'''
	backend = GrfSequence()
	golden_ast_plot(backend).grid()
	labels = [command["text"] for command in backend.take_commands()
	          if command["op"] == GrfOp.TEXT.value]
	assert labels
	assert not [label for label in labels if "%" in label]


def test_grid_right_ascension_labels_are_in_hours():
	'''
	AST's own default RA format is hours (``2:00``). v1 set ``Format(1)=dms``,
	which actually forces *degrees* — a misunderstanding SPEC-09A does not
	carry forward, so v2 sets no Format attribute at all. This pins that.
	'''
	backend = GrfSequence()
	golden_ast_plot(backend).grid()
	labels = [command["text"] for command in backend.take_commands()
	          if command["op"] == GrfOp.TEXT.value]
	assert "2:00" in labels, f"expected an hours-format RA label among {labels}"


def test_recorded_grid_stream_is_a_valid_envelope_payload():
	from cornish.plot.envelope import LayerAction, make_envelope
	backend = GrfSequence()
	golden_ast_plot(backend).grid()
	make_envelope(layer_id="grid", action=LayerAction.REPLACE, kind="grid", z=-10,
	              commands=backend.take_commands(),
	              canvas={"gbox": [0.0, 0.0, 1.0, 1.0], "aspect": 1.0,
	                      "margin": 0.15, "char_height": 0.02})
