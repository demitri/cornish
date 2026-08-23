
'''
Run the JS conformance suite from the Python suite (SPEC-09A §7.4).

``cornish-view.js`` is half of a two-implementation contract: it must produce
the same geometry from the same fixtures as the Python writer, and reproduce
the same RFC 8785 canonical hashes. A JS suite that only runs when someone
remembers to type ``node --test`` is a convention, and conventions decay — so
the Python suite drives it, and a JS regression fails the build that a Python
developer actually runs.

Skipped (not silently passed) when node is unavailable.
'''

import os
import shutil
import subprocess

import pytest

STATIC_DIRECTORY = os.path.normpath(
	os.path.join(os.path.dirname(__file__), "..", "cornish", "view", "static"))


def test_js_conformance_suite_passes():
	node = shutil.which("node")
	if node is None:
		pytest.skip("node is not installed; the JS conformance suite cannot run here")

	result = subprocess.run([node, "--test"], cwd=STATIC_DIRECTORY,
	                        capture_output=True, text=True)
	if result.returncode != 0:
		pytest.fail("the cornish-view.js conformance suite failed:\n"
		            + result.stdout[-4000:] + "\n" + result.stderr[-2000:])

	# guard against a "pass" that ran nothing at all
	assert "# pass " in result.stdout, result.stdout[-2000:]
	passed = int(result.stdout.split("# pass ", 1)[1].split("\n", 1)[0])
	assert passed >= 30, f"only {passed} JS tests ran; the suite should cover far more"


def test_the_js_renderer_and_its_suite_are_present():
	for name in ("cornish-view.js", "cornish-view.test.mjs", "viewer.html",
	             "envelope.schema.json"):
		assert os.path.exists(os.path.join(STATIC_DIRECTORY, name)), f"missing static asset {name}"
