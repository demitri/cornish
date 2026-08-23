
'''
Colour handling for the plotting stack (SPEC-09A §1.1, §2.6).

Two deliberately distinct grammars live here:

* :func:`normalize_colour` — for **user colour parameters** (a draw method's
  ``colour=``, an ``attr`` command's colour string). Accepts a base-palette
  index, ``#rgb``/``#rrggbb``, or a CSS named colour.
* :func:`normalize_palette_entry` — for **palette entries** themselves
  (``canvas.palette``, ``SkyPlot(palette=...)``). Wider: it additionally
  accepts ``transparent`` and functional ``rgb()``/``rgba()`` forms, because
  palette slot 0 is legitimately transparent. It does *not* accept integers —
  a palette entry cannot reference the palette.

Neither ever falls back to a default: anything outside its grammar raises
``ValueError``. :data:`DEFAULT_PALETTE` passes :func:`normalize_palette_entry`
by construction and deliberately does *not* pass :func:`normalize_colour`.
'''

import re
from typing import Mapping, Sequence, Union

__all__ = ['DEFAULT_PALETTE', 'CSS_NAMED_COLOURS', 'normalize_colour', 'normalize_palette_entry']

#: The frozen 16-entry base palette (SPEC-09A §2.6, D20) — astplot-js
#: ``lightTheme.colours`` verbatim, so a base-palette index means the same
#: colour in every consumer of an envelope and theme switching works for free.
#: Index 0 is the background slot (transparent), index 1 the foreground;
#: 2–15 are the matplotlib "tab" cycle astplot-js already chose.
DEFAULT_PALETTE = (
	"rgba(0,0,0,0)", "#000000", "#d62728", "#2ca02c",
	"#1f77b4", "#17becf", "#9467bd", "#bcbd22",
	"#ff7f0e", "#8c564b", "#e377c2", "#7f7f7f",
	"#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5",
)

#: CSS Color Module Level 4 named colours -> ``#rrggbb`` (148 entries).
CSS_NAMED_COLOURS: Mapping[str, str] = {
	"aliceblue": "#f0f8ff", "antiquewhite": "#faebd7", "aqua": "#00ffff", "aquamarine": "#7fffd4",
	"azure": "#f0ffff", "beige": "#f5f5dc", "bisque": "#ffe4c4", "black": "#000000",
	"blanchedalmond": "#ffebcd", "blue": "#0000ff", "blueviolet": "#8a2be2", "brown": "#a52a2a",
	"burlywood": "#deb887", "cadetblue": "#5f9ea0", "chartreuse": "#7fff00",
	"chocolate": "#d2691e", "coral": "#ff7f50", "cornflowerblue": "#6495ed",
	"cornsilk": "#fff8dc", "crimson": "#dc143c", "cyan": "#00ffff", "darkblue": "#00008b",
	"darkcyan": "#008b8b", "darkgoldenrod": "#b8860b", "darkgray": "#a9a9a9",
	"darkgreen": "#006400", "darkgrey": "#a9a9a9", "darkkhaki": "#bdb76b",
	"darkmagenta": "#8b008b", "darkolivegreen": "#556b2f", "darkorange": "#ff8c00",
	"darkorchid": "#9932cc", "darkred": "#8b0000", "darksalmon": "#e9967a",
	"darkseagreen": "#8fbc8f", "darkslateblue": "#483d8b", "darkslategray": "#2f4f4f",
	"darkslategrey": "#2f4f4f", "darkturquoise": "#00ced1", "darkviolet": "#9400d3",
	"deeppink": "#ff1493", "deepskyblue": "#00bfff", "dimgray": "#696969", "dimgrey": "#696969",
	"dodgerblue": "#1e90ff", "firebrick": "#b22222", "floralwhite": "#fffaf0",
	"forestgreen": "#228b22", "fuchsia": "#ff00ff", "gainsboro": "#dcdcdc",
	"ghostwhite": "#f8f8ff", "gold": "#ffd700", "goldenrod": "#daa520", "gray": "#808080",
	"green": "#008000", "greenyellow": "#adff2f", "grey": "#808080", "honeydew": "#f0fff0",
	"hotpink": "#ff69b4", "indianred": "#cd5c5c", "indigo": "#4b0082", "ivory": "#fffff0",
	"khaki": "#f0e68c", "lavender": "#e6e6fa", "lavenderblush": "#fff0f5", "lawngreen": "#7cfc00",
	"lemonchiffon": "#fffacd", "lightblue": "#add8e6", "lightcoral": "#f08080",
	"lightcyan": "#e0ffff", "lightgoldenrodyellow": "#fafad2", "lightgray": "#d3d3d3",
	"lightgreen": "#90ee90", "lightgrey": "#d3d3d3", "lightpink": "#ffb6c1",
	"lightsalmon": "#ffa07a", "lightseagreen": "#20b2aa", "lightskyblue": "#87cefa",
	"lightslategray": "#778899", "lightslategrey": "#778899", "lightsteelblue": "#b0c4de",
	"lightyellow": "#ffffe0", "lime": "#00ff00", "limegreen": "#32cd32", "linen": "#faf0e6",
	"magenta": "#ff00ff", "maroon": "#800000", "mediumaquamarine": "#66cdaa",
	"mediumblue": "#0000cd", "mediumorchid": "#ba55d3", "mediumpurple": "#9370db",
	"mediumseagreen": "#3cb371", "mediumslateblue": "#7b68ee", "mediumspringgreen": "#00fa9a",
	"mediumturquoise": "#48d1cc", "mediumvioletred": "#c71585", "midnightblue": "#191970",
	"mintcream": "#f5fffa", "mistyrose": "#ffe4e1", "moccasin": "#ffe4b5",
	"navajowhite": "#ffdead", "navy": "#000080", "oldlace": "#fdf5e6", "olive": "#808000",
	"olivedrab": "#6b8e23", "orange": "#ffa500", "orangered": "#ff4500", "orchid": "#da70d6",
	"palegoldenrod": "#eee8aa", "palegreen": "#98fb98", "paleturquoise": "#afeeee",
	"palevioletred": "#db7093", "papayawhip": "#ffefd5", "peachpuff": "#ffdab9",
	"peru": "#cd853f", "pink": "#ffc0cb", "plum": "#dda0dd", "powderblue": "#b0e0e6",
	"purple": "#800080", "rebeccapurple": "#663399", "red": "#ff0000", "rosybrown": "#bc8f8f",
	"royalblue": "#4169e1", "saddlebrown": "#8b4513", "salmon": "#fa8072",
	"sandybrown": "#f4a460", "seagreen": "#2e8b57", "seashell": "#fff5ee", "sienna": "#a0522d",
	"silver": "#c0c0c0", "skyblue": "#87ceeb", "slateblue": "#6a5acd", "slategray": "#708090",
	"slategrey": "#708090", "snow": "#fffafa", "springgreen": "#00ff7f", "steelblue": "#4682b4",
	"tan": "#d2b48c", "teal": "#008080", "thistle": "#d8bfd8", "tomato": "#ff6347",
	"turquoise": "#40e0d0", "violet": "#ee82ee", "wheat": "#f5deb3", "white": "#ffffff",
	"whitesmoke": "#f5f5f5", "yellow": "#ffff00", "yellowgreen": "#9acd32",
}

_HEX3 = re.compile(r"^#[0-9a-f]{3}$")
_HEX6 = re.compile(r"^#[0-9a-f]{6}$")
# rgb()/rgba() functional syntax: 3 or 4 comma- or space-separated components,
# each an integer, a percentage, or a decimal alpha. This is a SYNTAX check
# only — the consumer of a palette entry is a CSS engine, which owns the
# semantics (SPEC-09A §1.1: "passed through verbatim after a syntax check").
_RGB_FUNC = re.compile(
	r"^rgba?\(\s*"
	r"(-?\d+(?:\.\d+)?%?)\s*[,\s]\s*"
	r"(-?\d+(?:\.\d+)?%?)\s*[,\s]\s*"
	r"(-?\d+(?:\.\d+)?%?)"
	r"(?:\s*[,/]\s*|\s+)?(-?\d+(?:\.\d+)?%?)?"
	r"\s*\)$")


def _hex_or_none(text: str):
	'''
	Return the lower-cased, 6-digit form of a ``#rgb``/``#rrggbb`` string, or
	``None`` when *text* is not a hex colour at all.

	``None`` means "not this form", never "use a default" — every caller turns
	an exhausted set of forms into a ``ValueError``.
	'''
	lowered = text.strip().lower()
	if _HEX6.match(lowered):
		return lowered
	if _HEX3.match(lowered):
		return "#" + "".join(character * 2 for character in lowered[1:])
	return None


def normalize_colour(value: Union[str, int], palette: Sequence[str]) -> Union[str, int]:
	'''
	Normalize a USER colour parameter (SPEC-09A §1.1).

	:param value: a base-palette index, ``#rgb``/``#rrggbb``, or a CSS named colour
	:param palette: the palette an integer index is validated against
	:returns: the validated ``int`` index unchanged, or a lower-cased ``#rrggbb`` string
	:raises TypeError: if *value* is neither an ``int`` nor a ``str``
	:raises ValueError: if *value* is an out-of-range index or an unrecognized colour
	'''
	if isinstance(value, bool):
		# bool is an int subclass; a True/False colour index is always a mistake
		raise TypeError(f"A colour must be a palette index (int) or a CSS colour string; got {value!r}.")
	if isinstance(value, int):
		if not 0 <= value < len(palette):
			raise ValueError(f"Colour index {value} is outside the palette (0..{len(palette) - 1}).")
		return value
	if not isinstance(value, str):
		raise TypeError("A colour must be a palette index (int) or a CSS colour string; "
		                f"got {type(value).__name__}.")

	hex_form = _hex_or_none(value)
	if hex_form is not None:
		return hex_form
	named = CSS_NAMED_COLOURS.get(value.strip().lower())
	if named is not None:
		return named
	raise ValueError(f"'{value}' is not a recognized colour: expected a palette index, "
	                 "'#rgb'/'#rrggbb', or a CSS named colour.")


def normalize_palette_entry(value: str) -> str:
	'''
	Normalize a PALETTE entry — the deliberately wider grammar (SPEC-09A §1.1).

	:param value: ``#rgb``/``#rrggbb``, a CSS named colour, ``transparent``, or an ``rgb()``/``rgba()`` form
	:returns: the normalized entry (hex lower-cased and expanded; functional forms case-folded and passed through)
	:raises TypeError: if *value* is not a ``str`` — an ``int`` is rejected outright, since a palette entry cannot reference the palette
	:raises ValueError: if *value* is not one of the accepted forms
	'''
	if not isinstance(value, str):
		raise TypeError("A palette entry must be a CSS colour string "
		                f"(an int cannot reference the palette); got {type(value).__name__}.")
	text = value.strip()
	lowered = text.lower()
	if lowered == "transparent":
		return "transparent"
	hex_form = _hex_or_none(text)
	if hex_form is not None:
		return hex_form
	named = CSS_NAMED_COLOURS.get(lowered)
	if named is not None:
		return named
	if _RGB_FUNC.match(lowered):
		return lowered
	raise ValueError(f"'{value}' is not a valid palette entry: expected '#rgb'/'#rrggbb', "
	                 "a CSS named colour, 'transparent', or an rgb()/rgba() form.")
