
'''
Golden draw-stream fixtures (SPEC-09A §7.3, §7.6).

These fixtures are the **canonical** copy of the cross-implementation contract:
astplot-js vendors them, replays them through its own renderer, reproduces
their canonical hashes with its own ``canonicalJSON``, and — the
cross-implementation prize — drives its own AST-WASM Grf bridge with golden 1's
header and passes the same tolerant compare.

The compare is deliberately **tolerant**, because a golden produced by one AST
build is checked against another: label text must match exactly (that is AST's
tick *decision*, which must not drift), command counts within ±10% (polyline
point counts legitimately vary between AST versions), and every coordinate must
be finite and inside the viewport. Float equality is NOT asserted.

Regenerating requires ``--regenerate-goldens``; nothing here ever refreshes a
fixture silently.
'''

import json
import os
from collections import Counter

import astropy.units as u
import numpy as np
import pytest

from cornish import ASTCircle, ASTICRSFrame
from cornish.plot.envelope import canonical_json, envelope_hash, validate_document
from cornish.plot.svg import render_svg
from cornish.plot.sky import SkyPlot

FIXTURE_DIRECTORY = os.path.join(os.path.dirname(__file__), "fixtures", "draw_streams")

#: The tolerance on per-op command counts between AST builds.
COUNT_TOLERANCE = 0.10


def build_golden1() -> dict:
	''' The PoC TAN field (CRVAL 30/45, 2° radius), grid only. '''
	extent = ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg)
	return SkyPlot(extent=extent, plot_id="golden1").document()


def build_golden2() -> dict:
	''' TAN + a circle region outline + 3 markers — the §9 curl example. '''
	extent = ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=2 * u.deg)
	plot = SkyPlot(extent=extent, plot_id="golden2")
	plot.addRegionOutline(ASTCircle(frame=ASTICRSFrame(), center=[30.0, 45.0], radius=1 * u.deg))
	plot.addPoints(ra=[29.5, 30.0, 30.5], dec=[44.8, 45.0, 45.2])
	return plot.document()


def build_golden3() -> dict:
	''' All-sky AIT, galactic. '''
	return SkyPlot(all_sky=True, system="GALACTIC", plot_id="golden3").document()


GOLDEN_BUILDERS = {
	"golden1": build_golden1,
	"golden2": build_golden2,
	"golden3": build_golden3,
}


def fixture_path(name: str) -> str:
	return os.path.join(FIXTURE_DIRECTORY, f"{name}.json")


def hash_path(name: str) -> str:
	return os.path.join(FIXTURE_DIRECTORY, f"{name}.hash")


def load_golden(name: str) -> dict:
	with open(fixture_path(name), encoding="utf-8") as stream:
		return json.load(stream)


def summarize(document: dict) -> dict:
	''' The tolerant-compare summary of a document. '''
	operations = Counter()
	labels = Counter()
	layers = []
	for layer in document["layers"]:
		layers.append((layer["id"], layer.get("kind"), layer.get("z", 0)))
		for command in layer["commands"]:
			operations[command["op"]] += 1
			if command["op"] == "text":
				labels[command["text"]] += 1
	return {"operations": operations, "labels": labels, "layers": layers}


def viewport(document: dict):
	canvas = document["canvas"]
	x1, y1, x2, y2 = (float(value) for value in canvas["gbox"])
	bleed = float(canvas.get("margin", 0.15)) * max(x2 - x1, y2 - y1)
	return x1 - bleed, y1 - bleed, x2 + bleed, y2 + bleed


# ---------------------------------------------------------------------------
# regeneration
# ---------------------------------------------------------------------------

def test_regenerate_goldens(request):
	'''
	Rewrites the fixtures — but only under ``--regenerate-goldens``, and it
	prints a diff summary of what changed rather than refreshing silently.
	'''
	if not request.config.getoption("--regenerate-goldens"):
		pytest.skip("pass --regenerate-goldens to rewrite the golden fixtures")

	os.makedirs(FIXTURE_DIRECTORY, exist_ok=True)
	for name, builder in GOLDEN_BUILDERS.items():
		document = builder()
		validate_document(document)
		serialized = canonical_json(document)
		digest = envelope_hash(document)

		if os.path.exists(fixture_path(name)):
			previous = load_golden(name)
			before, after = summarize(previous), summarize(document)
			print(f"\n[{name}] ops {dict(before['operations'])} -> {dict(after['operations'])}")
			added = set(after["labels"]) - set(before["labels"])
			removed = set(before["labels"]) - set(after["labels"])
			print(f"[{name}] labels added {sorted(added)} removed {sorted(removed)}")
			with open(hash_path(name), encoding="utf-8") as stream:
				print(f"[{name}] hash {stream.read().strip()[:16]} -> {digest[:16]}")
		else:
			print(f"\n[{name}] NEW fixture, hash {digest[:16]}")

		with open(fixture_path(name), "w", encoding="utf-8") as stream:
			stream.write(serialized)
		with open(hash_path(name), "w", encoding="utf-8") as stream:
			stream.write(digest + "\n")


# ---------------------------------------------------------------------------
# the fixtures themselves
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_fixture_exists(name):
	assert os.path.exists(fixture_path(name)), \
		f"missing golden fixture {name}.json (regenerate with --regenerate-goldens)"
	assert os.path.exists(hash_path(name)), f"missing golden hash {name}.hash"


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_fixture_is_a_valid_document(name):
	validate_document(load_golden(name))


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_fixture_hash_is_pinned(name):
	'''
	The cross-language pin: astplot-js reproduces these with its own
	``canonicalJSON``. A change here is a contract change, not a refactor.
	'''
	with open(hash_path(name), encoding="utf-8") as stream:
		expected = stream.read().strip()
	assert envelope_hash(load_golden(name)) == expected


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_fixture_file_is_canonical_json(name):
	''' The stored bytes are themselves the canonical form. '''
	with open(fixture_path(name), encoding="utf-8") as stream:
		stored = stream.read()
	assert stored == canonical_json(load_golden(name))


# ---------------------------------------------------------------------------
# tolerant compare against a live rebuild
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_live_rebuild_matches_the_golden(name):
	'''
	The tolerant compare (§7.3). Labels must match EXACTLY — they encode AST's
	tick decisions, and a drift there means the layout changed. Counts are
	allowed ±10% because polyline point counts legitimately vary across AST
	versions.
	'''
	stored = summarize(load_golden(name))
	live = summarize(GOLDEN_BUILDERS[name]())

	assert live["labels"] == stored["labels"], (
		f"{name}: label set changed\n  only live: "
		f"{sorted(set(live['labels']) - set(stored['labels']))}\n  only stored: "
		f"{sorted(set(stored['labels']) - set(live['labels']))}")

	assert live["layers"] == stored["layers"], f"{name}: layer identity/order changed"

	for operation, stored_count in stored["operations"].items():
		live_count = live["operations"][operation]
		allowed = max(1, COUNT_TOLERANCE * stored_count)
		assert abs(live_count - stored_count) <= allowed, (
			f"{name}: '{operation}' count {live_count} is outside "
			f"{stored_count} +/- {allowed:.1f}")


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_every_coordinate_is_finite_and_inside_the_viewport(name):
	document = load_golden(name)
	left, bottom, right, top = viewport(document)
	for layer in document["layers"]:
		for command in layer["commands"]:
			if command["op"] in ("line", "mark"):
				points = zip(command["x"], command["y"])
			elif command["op"] == "text":
				points = [(command["x"], command["y"])]
			else:
				continue
			for x, y in points:
				assert np.isfinite(x) and np.isfinite(y), f"{name}/{layer['id']}: non-finite point"
				assert left <= x <= right and bottom <= y <= top, (
					f"{name}/{layer['id']}: ({x}, {y}) is outside the viewport "
					f"[{left}, {right}] x [{bottom}, {top}]")


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_every_golden_renders(name):
	''' §7.4: every golden renders, parses, and contains no NaN. '''
	import xml.etree.ElementTree as ElementTree
	svg = render_svg(load_golden(name))
	ElementTree.fromstring(svg)
	assert "NaN" not in svg


@pytest.mark.parametrize("name", sorted(GOLDEN_BUILDERS))
def test_render_is_byte_stable_for_every_golden(name):
	document = load_golden(name)
	assert render_svg(document) == render_svg(document)


def test_golden1_labels_are_hours_for_right_ascension():
	''' AST's default RA format; v2 sets no Format attribute (§2.2). '''
	labels = summarize(load_golden("golden1"))["labels"]
	assert "2:00" in labels


def test_golden2_has_a_region_and_a_points_layer():
	layers = {layer["id"]: layer.get("kind") for layer in load_golden("golden2")["layers"]}
	assert "region-1" in layers and layers["region-1"] == "region"
	assert "points-2" in layers and layers["points-2"] == "points"


def test_golden2_marks_exactly_three_points():
	document = load_golden("golden2")
	points = [layer for layer in document["layers"] if layer["id"] == "points-2"][0]
	marks = [command for command in points["commands"] if command["op"] == "mark"]
	assert len(marks) == 1 and len(marks[0]["x"]) == 3


def test_golden3_is_all_sky_galactic():
	document = load_golden("golden3")
	assert document["canvas"]["aspect"] == 0.5
	labels = summarize(document)["labels"]
	assert any("Galactic" in label or "galactic" in label for label in labels), \
		f"expected a galactic title among {sorted(labels)}"


def test_fixtures_readme_exists():
	''' §7.6: the fixture directory documents its own compare and regeneration rules. '''
	assert os.path.exists(os.path.join(FIXTURE_DIRECTORY, "FIXTURES.md"))
