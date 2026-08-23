#!/usr/bin/env python3
'''
Dev-time generator for the shipped font-metrics table (SPEC-09A §6.1).

The table lets ``cornish.plot.metrics.FontMetrics`` answer AST's ``Qch``/
``TxExt`` questions without a renderer, which is what makes SkyPlot v2's
layout identical in an SVG file, a notebook, a browser, and a matplotlib
window (D17). Neither fontTools nor matplotlib is a runtime dependency of
cornish — the generated JSON ships in the wheel; this script runs only when
the table format or the coverage set changes.

    python3 generate_font_metrics.py [--font /path/to/DejaVuSans.ttf] [-o OUT]

With no ``--font``, the script falls back to matplotlib's bundled copy.

**Fails loud on any glyph in the coverage set that the font does not
provide** — a silently-missing glyph would degrade layout accuracy
invisibly, which is exactly the failure the coverage set exists to prevent.
'''

import argparse
import json
import os
import sys


#: The frozen coverage set (SPEC-09A §6.2). AST's own axis labels are plain
#: ASCII [verified], but user ``addText`` strings in astronomy are Greek-heavy
#: (α, δ, λ, M⊙), and every uncovered glyph costs layout accuracy plus a
#: per-glyph diagnostic at render time.
def coverage_set() -> list:
	characters = [chr(code) for code in range(0x20, 0x7F)]          # U+0020–U+007E
	characters += ["°", "′", "″"]                     # ° ′ ″
	characters += [chr(code) for code in range(0x391, 0x3AA) if code != 0x3A2]  # Α–Ω (U+03A2 reserved)
	characters += [chr(code) for code in range(0x3B1, 0x3CA)]        # α–ω
	characters += ["Å", "µ", "±", "×", "·", "⊙", "–", "—"]
	return characters


def default_font_path() -> str:
	'''
	Locate matplotlib's bundled DejaVuSans.ttf. Raises rather than returning a
	sentinel: a caller with no ``--font`` and no matplotlib has no font at all.
	'''
	try:
		import matplotlib
	except ImportError as error:
		raise SystemExit("No --font given and matplotlib is not importable, so its bundled "
		                 "DejaVuSans.ttf cannot be located. Pass --font PATH.") from error
	path = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans.ttf")
	if not os.path.exists(path):
		raise SystemExit(f"matplotlib is installed but {path} does not exist. Pass --font PATH.")
	return path


def build_table(font_path: str) -> dict:
	'''
	Read *font_path* with fontTools and return the §6.2 table dict.

	:raises SystemExit: if any glyph of the coverage set is absent from the font
	'''
	try:
		from fontTools.ttLib import TTFont
	except ImportError as error:
		raise SystemExit("fontTools is required to generate the metrics table "
		                 "(it is NOT a cornish runtime dependency): pip install fonttools") from error

	font = TTFont(font_path)
	units_per_em = font["head"].unitsPerEm
	cmap = font.getBestCmap()
	hmtx = font["hmtx"]

	missing = []
	advance = {}
	for character in coverage_set():
		glyph_name = cmap.get(ord(character))
		if glyph_name is None or glyph_name not in hmtx.metrics:
			missing.append(f"U+{ord(character):04X} ({character!r})")
			continue
		advance[character] = round(hmtx.metrics[glyph_name][0] / units_per_em, 4)

	if missing:
		# fail loud: a partial table would silently mis-measure exactly the
		# strings the coverage set was extended to cover
		raise SystemExit(f"{len(missing)} glyph(s) of the coverage set are missing from "
		                 f"{font_path}:\n  " + "\n  ".join(missing))

	if "0" not in advance:
		raise SystemExit("The coverage set must contain '0' — its advance is the default_advance.")

	return {
		"schema": 1,
		"font_family": "DejaVu Sans",
		"ascent": round(font["hhea"].ascent / units_per_em, 4),
		"descent": round(font["hhea"].descent / units_per_em, 4),
		"advance": advance,
		"default_advance": advance["0"],
	}


def main(argv=None) -> None:
	parser = argparse.ArgumentParser(description="Generate the cornish font-metrics table (SPEC-09A §6).")
	parser.add_argument("--font", default=None, help="path to DejaVuSans.ttf (default: matplotlib's bundled copy)")
	parser.add_argument("-o", "--output", default=None,
	                    help="output path (default: cornish/plot/data/dejavu_sans_metrics.json)")
	arguments = parser.parse_args(argv)

	font_path = arguments.font or default_font_path()
	table = build_table(font_path)

	output = arguments.output
	if output is None:
		here = os.path.dirname(os.path.abspath(__file__))
		output = os.path.join(here, "..", "cornish", "plot", "data", "dejavu_sans_metrics.json")
		output = os.path.normpath(output)
	os.makedirs(os.path.dirname(output), exist_ok=True)
	with open(output, "w", encoding="utf-8") as stream:
		json.dump(table, stream, ensure_ascii=False, indent=1, sort_keys=True)
		stream.write("\n")
	print(f"wrote {output}: {len(table['advance'])} glyphs from {font_path}", file=sys.stderr)


if __name__ == "__main__":
	main()
