
'''
``GrfSequence`` — the recording Grf backend (SPEC-09A §1.1, D17).

This is the **only** Grf backend cornish implements. ``Ast.Plot`` draws into it
exactly as it would into a live renderer; instead of putting ink on a canvas it
records the primitive stream, which every output path (SVG, JSON, browser,
matplotlib) then *replays*. One backend means one metrics source and one
geometry truth, so the same plot is identical everywhere — which is not true of
``grf_matplotlib``, whose ``TxExt`` measures by drawing into the live figure and
therefore lays labels out differently per backend and per DPI.

Coordinates are recorded exactly as AST supplies them, in gbox units.

Attribute chatter is coalesced. AST sets an attribute, draws, and restores it
constantly; emitting one command per ``Attr`` call would triple the stream for
no information. Instead the *state* is tracked and a net ``attr`` command is
emitted lazily, immediately before the next draw op, and only for the fields
whose effective value differs from what was last emitted.
'''

import logging
import math
from typing import Optional, Sequence

import starlink.Ast as Ast

from .envelope import ATTR_DEFAULTS, GrfOp, PrimType
from .metrics import FontMetrics
from .palette import DEFAULT_PALETTE, normalize_palette_entry

__all__ = ['GrfSequence']

logger = logging.getLogger("cornish")

#: AST's Grf attribute codes -> the envelope's ``attr`` field names (§3.3).
_ATTR_NAMES = {
	Ast.grfSTYLE: "style",
	Ast.grfWIDTH: "width",
	Ast.grfSIZE: "size",
	Ast.grfFONT: "font",
	Ast.grfCOLOUR: "colour",
}

#: AST's Grf primitive codes -> the envelope's ``prim`` values (§3.3).
_PRIM_NAMES = {
	Ast.grfLINE: PrimType.LINE.value,
	Ast.grfMARK: PrimType.MARK.value,
	Ast.grfTEXT: PrimType.TEXT.value,
}

#: Which attributes the schema gives meaning to for each primitive class
#: (§3.3: width applies to lines/marks, style to lines, size to marks/text,
#: font to text, colour to everything). AST is free to set any attribute on any
#: primitive and the full state is tracked faithfully for the query protocol;
#: this table governs only what gets SERIALIZED, so the stream carries no
#: fields the schema declares meaningless for that primitive. It is a
#: documented projection onto the schema's vocabulary, not an error path.
_PRIM_VOCABULARY = {
	Ast.grfLINE: ("colour", "width", "style"),
	Ast.grfMARK: ("colour", "width", "size"),
	Ast.grfTEXT: ("colour", "size", "font"),
}

#: Attributes that are integers in the envelope (the rest are float multipliers).
#: ``colour`` belongs here: AST hands it over as a double like every other Grf
#: attribute, but it is a palette *index*, and recording it as 11.0 would both
#: violate the schema and defeat the index/CSS-string resolution in
#: :meth:`GrfSequence.take_commands`.
_INTEGER_ATTRS = ("style", "font", "colour")


class GrfSequence:
	'''
	A recording Grf backend for :class:`starlink.Ast.Plot`.

	Implements the full Grf contract AST requires — ``Line``/``Mark``/``Text``/
	``TxExt``/``Qch``/``Scales``/``Attr``/``Cap``/``BBuf``/``EBuf``/``Flush`` —
	plus pyast's ``ColToInt``/``IntToCol``, which are **mandatory**: without them
	``plot.set("Colour(Grid)=red")`` raises ``BADAT``.

	:param char_height: height of size-1.0 text in gbox y units (the ``Qch`` answer)
	:param palette: base palette; ``None`` uses :data:`~cornish.plot.palette.DEFAULT_PALETTE`
	:param metrics: font metrics; ``None`` loads the shipped table
	:raises ValueError: if *char_height* is not positive and finite, or *palette* is malformed
	'''

	def __init__(self, *, char_height: float = 0.02,
	             palette: Optional[Sequence[str]] = None,
	             metrics: Optional[FontMetrics] = None) -> None:
		char_height = float(char_height)
		if not math.isfinite(char_height) or char_height <= 0:
			raise ValueError(f"char_height must be a positive finite number; got {char_height!r}.")
		self._char_height = char_height

		if palette is None:
			entries = list(DEFAULT_PALETTE)
		else:
			entries = [normalize_palette_entry(entry) for entry in palette]
			if len(entries) < 2:
				raise ValueError("A palette needs at least 2 entries "
				                 f"(0 = background, 1 = foreground); got {len(entries)}.")
		self._palette = entries
		#: indices below this are theme slots and serialize as ints; at or above
		#: it are colours grown by ColToInt, which serialize as literal CSS
		self._base_palette_size = len(entries)

		self._metrics = metrics if metrics is not None else FontMetrics()

		self._commands = []
		# full per-primitive attribute state, as AST believes it to be
		self._state = {prim: dict(ATTR_DEFAULTS) for prim in _PRIM_NAMES}
		# what the stream has actually said since the last take_commands()
		self._emitted = {prim: dict(ATTR_DEFAULTS) for prim in _PRIM_NAMES}

	# ------------------------------------------------------------------
	# properties
	# ------------------------------------------------------------------

	@property
	def palette(self) -> list:
		''' The palette, base entries plus any grown by ``ColToInt`` (a copy). '''
		return list(self._palette)

	@property
	def char_height(self) -> float:
		''' Height of size-1.0 text, in gbox y units. '''
		return self._char_height

	@property
	def metrics(self) -> FontMetrics:
		''' The font metrics answering ``Qch``/``TxExt``. '''
		return self._metrics

	# ------------------------------------------------------------------
	# recording internals
	# ------------------------------------------------------------------

	def _coordinates(self, values, context: str) -> list:
		'''
		Coerce an AST coordinate sequence to a list of floats, refusing
		``Ast.BAD`` and non-finite values.

		``Ast.BAD`` is a *finite* float (-1.797e308), so it would sail through
		any ordinary finiteness check and be recorded as a real coordinate —
		producing a stream that draws a line to the edge of the universe. AST
		breaks its own polylines at bad points and should never hand one over;
		if it ever does, that is a bug worth seeing, so this raises rather than
		dropping the point (dropping would silently deform the geometry).
		'''
		coordinates = []
		for index, value in enumerate(values):
			number = float(value)
			if number == Ast.BAD:
				raise ValueError(f"{context}: AST supplied Ast.BAD at index {index}. A bad "
				                 "coordinate must never enter a draw stream; AST is expected to "
				                 "break polylines at bad points itself.")
			if not math.isfinite(number):
				raise ValueError(f"{context}: coordinate {number!r} at index {index} is not finite.")
			coordinates.append(number)
		return coordinates

	def _flush_attributes(self, prim) -> None:
		'''
		Emit the net ``attr`` command for *prim*, if its effective state has
		drifted from what the stream last said. Called immediately before every
		draw op — never on ``Attr`` itself, which is where the chatter is.
		'''
		vocabulary = _PRIM_VOCABULARY[prim]
		current, emitted = self._state[prim], self._emitted[prim]
		changed = {name: current[name] for name in vocabulary if current[name] != emitted[name]}
		if not changed:
			return
		command = {"op": GrfOp.ATTR.value, "prim": _PRIM_NAMES[prim]}
		for name, value in changed.items():
			command[name] = int(value) if name in _INTEGER_ATTRS else float(value)
			emitted[name] = current[name]
		self._commands.append(command)

	# ------------------------------------------------------------------
	# the Grf contract (called by AST; not user API)
	# ------------------------------------------------------------------

	def Line(self, n, x, y) -> None:
		''' Draw an *n*-point polyline. '''
		if n < 2:
			return                      # AST's own contract: fewer than 2 points draws nothing
		self._flush_attributes(Ast.grfLINE)
		self._commands.append({
			"op": GrfOp.LINE.value,
			"x": self._coordinates(x[:n], "Line"),
			"y": self._coordinates(y[:n], "Line"),
		})

	def Mark(self, n, x, y, type) -> None:
		''' Draw *n* markers of AST marker code *type*. '''
		if n < 1:
			return
		self._flush_attributes(Ast.grfMARK)
		self._commands.append({
			"op": GrfOp.MARK.value,
			"x": self._coordinates(x[:n], "Mark"),
			"y": self._coordinates(y[:n], "Mark"),
			"marker": int(type),
		})

	def Text(self, text, x, y, just, upx, upy) -> None:
		''' Draw a text string anchored at (*x*, *y*). '''
		if text is None or text == "":
			return                      # nothing to draw; AST asks for empty labels routinely
		self._flush_attributes(Ast.grfTEXT)
		anchor_x, anchor_y = self._coordinates([x, y], "Text")
		self._commands.append({
			"op": GrfOp.TEXT.value,
			"text": str(text),
			"x": anchor_x,
			"y": anchor_y,
			"just": self._normalize_just(just),
			"up": [float(upx), float(upy)],
		})

	def TxExt(self, text, x, y, just, upx, upy) -> tuple:
		'''
		The bounding box of *text* as it would be drawn: 8 floats, x×4 then y×4.

		Answered from the shipped metrics table at the text primitive's current
		size — never by drawing.
		'''
		em = self._char_height * float(self._state[Ast.grfTEXT]["size"])
		return self._metrics.extents(str(text), x, y, self._normalize_just(just), upx, upy, em)

	def Qch(self) -> tuple:
		'''
		``(chv, chh)`` — the height of vertically- and horizontally-drawn
		characters. Equal here because the gbox is isotropic by construction
		(§2.2), which is also what makes :meth:`Scales` honest.
		'''
		return (self._char_height, self._char_height)

	def Scales(self) -> tuple:
		''' ``(1.0, 1.0)``: gbox units are isotropic by construction (§2.2). '''
		return (1.0, 1.0)

	def Attr(self, attr, value, prim) -> float:
		'''
		Set graphics attribute *attr* for primitive class *prim*, returning its
		previous value. ``value == Ast.BAD`` is AST's *query* protocol: report
		the current value without changing it.
		'''
		if attr not in _ATTR_NAMES:
			raise ValueError(f"Unknown Grf attribute code {attr!r}; "
			                 f"expected one of {sorted(_ATTR_NAMES)}.")
		if prim not in _PRIM_NAMES:
			raise ValueError(f"Unknown Grf primitive code {prim!r}; "
			                 f"expected one of {sorted(_PRIM_NAMES)}.")
		name = _ATTR_NAMES[attr]
		previous = self._state[prim][name]
		if value != Ast.BAD:
			self._state[prim][name] = int(value) if name in _INTEGER_ATTRS else float(value)
		return float(previous)

	def Cap(self, cap, value) -> int:
		'''
		Report a backend capability.

		``SCALES`` and ``MJUST`` are supported; ``ESC`` returns 0, which is what
		makes AST deliver label text **plain** rather than sprinkled with ``%``
		escape sequences — the envelope's text op depends on that (§3.3).
		'''
		if cap in (Ast.grfSCALES, Ast.grfMJUST):
			return 1
		if cap == Ast.grfESC:
			return 0
		return 0

	def BBuf(self) -> None:
		''' Begin buffering. A recorder has nothing to batch. '''
		return None

	def EBuf(self) -> None:
		''' End buffering. A recorder has nothing to flush. '''
		return None

	def Flush(self) -> None:
		''' Flush output. A recorder has nothing pending. '''
		return None

	def ColToInt(self, colour) -> int:
		'''
		Resolve a colour to a palette index, growing the palette for a colour
		it has not seen. An unrecognized colour returns −1, which is how the
		Grf contract says "not a legal colour" — AST then reports the error
		itself, so nothing is silently substituted.
		'''
		if isinstance(colour, bool):
			return -1
		if isinstance(colour, int):
			return colour                    # int passthrough (§1.1)
		try:
			return int(colour)               # AST may hand an index over as a string
		except (TypeError, ValueError):
			pass
		try:
			normalized = normalize_palette_entry(colour)
		except (TypeError, ValueError):
			return -1
		for index, entry in enumerate(self._palette):
			if entry == normalized:
				return index
		self._palette.append(normalized)
		return len(self._palette) - 1

	def IntToCol(self, colour) -> Optional[str]:
		''' The CSS colour at palette index *colour*, or ``None`` if out of range. '''
		try:
			index = int(colour)
		except (TypeError, ValueError):
			return None
		if 0 <= index < len(self._palette):
			return self._palette[index]
		return None

	# ------------------------------------------------------------------
	# capture API
	# ------------------------------------------------------------------

	def take_commands(self) -> list:
		'''
		Pop every command recorded since the last take.

		The attr-state baseline resets to the canvas defaults afterwards, so
		each taken chunk is self-contained: the next layer re-emits any
		non-default state at its head. That is what lets a consumer delete or
		reorder one layer without restyling its neighbours (§3.4, D18).

		Colour indices are resolved here (D18): an index inside the base palette
		stays an ``int``, so consumers can re-theme it; an index grown by
		``ColToInt`` is written out as its literal CSS string, so the layer
		does not depend on a palette the consumer never received.
		'''
		commands = self._commands
		self._commands = []
		self._emitted = {prim: dict(ATTR_DEFAULTS) for prim in _PRIM_NAMES}
		return [self._resolve_colours(command) for command in commands]

	def _resolve_colours(self, command: dict) -> dict:
		if command.get("op") != GrfOp.ATTR.value or "colour" not in command:
			return command
		index = command["colour"]
		if not isinstance(index, int) or isinstance(index, bool):
			return command
		if 0 <= index < self._base_palette_size:
			return command                       # a theme slot: keep it re-themeable
		resolved = dict(command)
		if 0 <= index < len(self._palette):
			resolved["colour"] = self._palette[index]
		# an index outside the palette entirely is left as an int on purpose:
		# §3.3 gives it a visible render-time fallback (foreground + diagnostic),
		# which is strictly better than inventing a colour for it here
		return resolved

	@staticmethod
	def _normalize_just(just) -> str:
		'''
		Normalize AST's justification code to the schema's 2-character form.
		AST sends legal codes; anything else is a real bug and raises rather
		than being defaulted to a guess.
		'''
		text = str(just).strip().upper()
		if len(text) != 2 or text[0] not in "TCBM" or text[1] not in "LCR":
			raise ValueError(f"AST supplied the justification {just!r}, which is not a legal "
			                 "2-character code (vertical T/C/B/M then horizontal L/C/R).")
		return text
