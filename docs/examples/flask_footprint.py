
'''
Flask **pattern 1**: render a sky plot server-side and serve it as an SVG
(SPEC-09A §7.5; SPEC-09's recommended production shape).

This file is the documentation. It is executable, it is imported by
``source/tests/test_flask_example.py``, and it is deliberately short, because
pattern 1 really is this small: ``render_svg`` is a **pure function** of a
draw-stream document, so serving a plot needs no daemon, no SSE, no background
worker and no shared state.

Run it::

    pip install flask
    python3 docs/examples/flask_footprint.py
    # then open http://127.0.0.1:5000/footprint.svg?ra=30&dec=45&radius=2

Why this shape, and not the blueprint:

* Pattern 1 has **no worker-model caveats** — no route holds a worker open, and
  there is no per-process registry to get out of step across gunicorn workers.
  ``cornish.view.flask_blueprint`` is for *live* plots (SSE, incremental
  updates) and carries the §5.6 caveats; read them before reaching for it.
* The ETag comes from ``envelope_hash``, the RFC 8785 canonical hash of the
  document. ``render_svg`` is deterministic, so equal documents produce equal
  bytes and therefore equal ETags — which is what makes conditional requests
  correct rather than merely plausible.
* Nothing here imports a transport library into the render path; the plotting
  stack has a machine-enforced gate against exactly that.
'''

import astropy.units as u
from flask import Flask, Response, request

from cornish import ASTCircle, ASTICRSFrame
from cornish.plot import SkyPlot, envelope_hash, render_svg

app = Flask(__name__)

#: Bounds for the user-supplied radius. A request outside them is refused
#: rather than clamped: a silently-resized field is a plot that lies.
MINIMUM_RADIUS_DEG = 0.001
MAXIMUM_RADIUS_DEG = 90.0


def build_document(ra_deg: float, dec_deg: float, radius_deg: float,
                   width: int = 640) -> dict:
	'''
	Build the draw-stream document for one footprint.

	Stage 1 builds the ``SkyPlot`` directly. At stage 2 this same function
	becomes ``plot_spec_to_envelopes(spec)`` — the declarative endpoint — and
	everything downstream of it is unchanged, which is the point of the
	envelope contract.

	**Note the explicit ``plot_id``, and copy it.** Left to itself a ``SkyPlot``
	invents ``plot-<pid>-<n>`` from a process-wide counter, so two identical
	requests would produce documents differing in that one field — identical
	SVG bytes, but a different ``envelope_hash`` every time, and therefore an
	ETag that never matches and a cache that never hits. Any hash-based caching
	over ``document()`` must pin the plot id.
	'''
	extent = ASTCircle(frame=ASTICRSFrame(),
	                   center=[ra_deg, dec_deg],
	                   radius=radius_deg * u.deg)
	plot = SkyPlot(extent=extent, width=width, plot_id="footprint")
	plot.addRegionOutline(extent, name="footprint")
	return plot.document()


@app.route("/footprint.svg")
def footprint():
	'''
	``GET /footprint.svg?ra=&dec=&radius=&width=`` -> an inline-able SVG.

	Bad parameters produce 400 with a message naming the offender — never a
	default substituted behind the caller's back.
	'''
	try:
		ra_deg = float(request.args.get("ra", 0.0))
		dec_deg = float(request.args.get("dec", 0.0))
		radius_deg = float(request.args.get("radius", 1.0))
		width = int(request.args.get("width", 640))
	except ValueError as error:
		return Response(f"bad parameter: {error}", status=400, mimetype="text/plain")

	if not MINIMUM_RADIUS_DEG <= radius_deg <= MAXIMUM_RADIUS_DEG:
		return Response(f"radius must lie in [{MINIMUM_RADIUS_DEG}, {MAXIMUM_RADIUS_DEG}] "
		                f"degrees; got {radius_deg}", status=400, mimetype="text/plain")
	if not 64 <= width <= 4096:
		return Response(f"width must lie in [64, 4096]; got {width}",
		                status=400, mimetype="text/plain")

	try:
		document = build_document(ra_deg, dec_deg, radius_deg, width=width)
	except ValueError as error:
		# cornish refuses an ill-defined field loudly (a non-finite bounding
		# circle, an unsupported system); that is a client error, not a 500
		return Response(str(error), status=400, mimetype="text/plain")

	# width changes the bytes, so it must change the ETag
	etag = f'"{envelope_hash(document)}-w{width}"'
	if request.headers.get("If-None-Match") == etag:
		return Response(status=304, headers={"ETag": etag})

	svg = render_svg(document, width=width)
	return Response(svg, mimetype="image/svg+xml", headers={"ETag": etag})


@app.route("/")
def index():
	''' A page embedding the SVG, which is the whole point of pattern 1. '''
	return Response(
		"<!doctype html><html><head><meta charset='utf-8'>"
		"<title>cornish footprint</title></head><body>"
		"<h1>M31 field</h1>"
		"<img src='/footprint.svg?ra=10.68&dec=41.27&radius=1.5' alt='footprint' width='640'>"
		"</body></html>",
		mimetype="text/html")


if __name__ == "__main__":
	app.run(debug=False)
