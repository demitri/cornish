/**
 * JS conformance suite for cornish-view.js (SPEC-09A §7.4).
 *
 *     node --test cornish-view.test.mjs
 *
 * Zero dependencies, matching astplot-js's node-test convention, under a
 * minimal DOM shim. This is the JS half of the two-implementation contract:
 * the same fixtures, the same canonical hashes, the same frozen fallbacks.
 */

import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { readFileSync, readdirSync } from "node:fs";
import { dirname, join } from "node:path";
import test from "node:test";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const FIXTURE_DIRECTORY = join(HERE, "..", "..", "..", "tests", "fixtures", "draw_streams");

// ---------------------------------------------------------------------------
// minimal DOM shim
// ---------------------------------------------------------------------------

class ShimElement {
	constructor(name, ownerDocument) {
		this.nodeName = name;
		this.tagName = name;
		this.nodeType = 1;
		this.ownerDocument = ownerDocument;
		this.attributes = {};
		this.childNodes = [];
		this.parentNode = null;
		this._text = "";
	}
	get firstChild() { return this.childNodes[0] || null; }
	get textContent() { return this._text; }
	set textContent(value) { this._text = String(value); }
	setAttribute(name, value) { this.attributes[name] = String(value); }
	getAttribute(name) {
		return Object.prototype.hasOwnProperty.call(this.attributes, name) ? this.attributes[name] : null;
	}
	appendChild(child) {
		child.parentNode = this;
		this.childNodes.push(child);
		return child;
	}
	removeChild(child) {
		const index = this.childNodes.indexOf(child);
		if (index !== -1) this.childNodes.splice(index, 1);
		child.parentNode = null;
		return child;
	}
	insertBefore(child, reference) {
		const index = this.childNodes.indexOf(reference);
		child.parentNode = this;
		if (index === -1) this.childNodes.push(child);
		else this.childNodes.splice(index, 0, child);
		return child;
	}
	/** Enough of a selector engine for [data-layer-id="..."]. */
	querySelectorAll(selector) {
		const match = /^\[([a-z-]+)="(.*)"\]$/.exec(selector);
		const results = [];
		const walk = (node) => {
			if (match && node.getAttribute(match[1]) === match[2]) results.push(node);
			node.childNodes.forEach(walk);
		};
		this.childNodes.forEach(walk);
		return results;
	}
	descendants() {
		const results = [];
		const walk = (node) => { results.push(node); node.childNodes.forEach(walk); };
		this.childNodes.forEach(walk);
		return results;
	}
}

const shimDocument = {
	createElementNS(_namespace, name) { return new ShimElement(name, shimDocument); },
};

function makeSVG() {
	return new ShimElement("svg", shimDocument);
}

/** Capture console.warn for the duration of `body`. */
function captureWarnings(body) {
	const original = console.warn;
	const messages = [];
	console.warn = (...args) => messages.push(args.join(" "));
	try {
		body();
	} finally {
		console.warn = original;
	}
	return messages;
}

const { CornishView, EnvelopeFormatError, canonicalJSON } = await import("./cornish-view.js");

const CANVAS = { gbox: [0, 0, 1, 1], aspect: 1.0, margin: 0.15, char_height: 0.02, width: 640 };

function envelope(overrides = {}) {
	return Object.assign({
		cornish_view: "1.0",
		type: "draw_stream",
		canvas: Object.assign({}, CANVAS),
		layer: { id: "grid", action: "replace" },
		commands: [],
	}, overrides);
}

// ---------------------------------------------------------------------------
// canonicalJSON
// ---------------------------------------------------------------------------

test("canonicalJSON sorts keys and omits whitespace", () => {
	assert.equal(canonicalJSON({ b: 1, a: 2, A: 3 }), '{"A":3,"a":2,"b":1}');
});

test("canonicalJSON uses native ES number formatting", () => {
	// exactly the cases where Python's repr disagrees, which is why the Python
	// side implements Number::toString rather than using repr
	assert.equal(canonicalJSON({ a: 1.0 }), '{"a":1}');
	assert.equal(canonicalJSON({ a: 1e-7 }), '{"a":1e-7}');
	assert.equal(canonicalJSON({ a: 0.00001 }), '{"a":0.00001}');
	assert.equal(canonicalJSON({ a: -0.0 }), '{"a":0}');
	assert.equal(canonicalJSON({ a: 1e21 }), '{"a":1e+21}');
	assert.equal(canonicalJSON({ a: 1e20 }), '{"a":100000000000000000000}');
});

test("canonicalJSON rejects non-finite numbers", () => {
	assert.throws(() => canonicalJSON({ a: NaN }), EnvelopeFormatError);
	assert.throws(() => canonicalJSON({ a: Infinity }), EnvelopeFormatError);
});

test("canonicalJSON passes non-ASCII through", () => {
	assert.equal(canonicalJSON({ t: "α δ ⊙" }), '{"t":"α δ ⊙"}');
});

// ---------------------------------------------------------------------------
// the golden fixtures — the cross-language contract
// ---------------------------------------------------------------------------

const goldenNames = readdirSync(FIXTURE_DIRECTORY)
	.filter((name) => name.endsWith(".json"))
	.map((name) => name.replace(/\.json$/, ""))
	.sort();

test("the golden fixtures are present", () => {
	assert.ok(goldenNames.length >= 3, `expected at least 3 goldens, found ${goldenNames.length}`);
});

for (const name of goldenNames) {
	test(`golden ${name}: canonical hash reproduces the Python pin`, () => {
		const raw = readFileSync(join(FIXTURE_DIRECTORY, `${name}.json`), "utf8");
		const expected = readFileSync(join(FIXTURE_DIRECTORY, `${name}.hash`), "utf8").trim();
		const digest = createHash("sha256").update(canonicalJSON(JSON.parse(raw)), "utf8").digest("hex");
		assert.equal(digest, expected,
			`${name}: JS canonicalJSON disagrees with the Python pin — the RFC 8785 number grammar is the contract`);
	});

	test(`golden ${name}: the stored file is already canonical`, () => {
		const raw = readFileSync(join(FIXTURE_DIRECTORY, `${name}.json`), "utf8");
		assert.equal(raw, canonicalJSON(JSON.parse(raw)));
	});

	test(`golden ${name}: replays without error`, () => {
		const svg = makeSVG();
		const view = new CornishView(svg);
		const document_ = JSON.parse(readFileSync(join(FIXTURE_DIRECTORY, `${name}.json`), "utf8"));
		view.apply(document_);

		const groups = svg.descendants().filter((node) => node.nodeName === "g" && node.getAttribute("data-layer-id"));
		assert.equal(groups.length, document_.layers.length, `${name}: layer group count`);
		// §3.5: no id attribute on layer groups
		for (const group of groups) assert.equal(group.getAttribute("id"), null);
		assert.ok(view.toSVGString().startsWith("<svg"));
	});
}

// ---------------------------------------------------------------------------
// shape dispatch (§3.1) — replicated independently of the Python validators
// ---------------------------------------------------------------------------

test("apply accepts an envelope", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope());
	assert.deepEqual(Array.from(view.layers.keys()), ["grid"]);
});

test("apply accepts a document", () => {
	const view = new CornishView(makeSVG());
	view.apply({
		cornish_view: "1.0", type: "draw_stream", canvas: CANVAS,
		layers: [{ id: "grid", kind: "grid", z: -10, commands: [] }],
	});
	assert.deepEqual(Array.from(view.layers.keys()), ["grid"]);
});

test("apply rejects a hybrid body carrying both layer and layers", () => {
	const view = new CornishView(makeSVG());
	const hybrid = envelope({ layers: [] });
	assert.throws(() => view.apply(hybrid), EnvelopeFormatError);
});

test("apply accepts the canvas-only envelope", () => {
	const view = new CornishView(makeSVG());
	const body = envelope();
	delete body.layer;
	delete body.commands;
	view.apply(body);
	assert.equal(view.layers.size, 0);
});

test("apply rejects a canvas-only envelope carrying commands", () => {
	const view = new CornishView(makeSVG());
	const body = envelope();
	delete body.layer;
	assert.throws(() => view.apply(body), EnvelopeFormatError);
});

test("apply rejects an unknown major version and accepts an unknown minor", () => {
	const view = new CornishView(makeSVG());
	assert.throws(() => view.apply(envelope({ cornish_view: "2.0" })), EnvelopeFormatError);
	view.apply(envelope({ cornish_view: "1.9" }));   // §3.7: unknown minor is fine
});

test("apply rejects a wrong document type", () => {
	const view = new CornishView(makeSVG());
	assert.throws(() => view.apply(envelope({ type: "plot_spec" })), EnvelopeFormatError);
});

// ---------------------------------------------------------------------------
// the layer state machine (§3.4)
// ---------------------------------------------------------------------------

test("append to a missing layer is an error, never an auto-create", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope());
	assert.throws(
		() => view.apply(envelope({ layer: { id: "ghost", action: "append" }, commands: [] })),
		EnvelopeFormatError);
});

test("delete of a missing layer is an error, never a silent no-op", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope());
	const body = envelope({ layer: { id: "ghost", action: "delete" } });
	delete body.commands;
	assert.throws(() => view.apply(body), EnvelopeFormatError);
});

test("append extends the stored commands", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope({ commands: [{ op: "line", x: [0, 1], y: [0, 1] }] }));
	view.apply(envelope({
		layer: { id: "grid", action: "append" },
		commands: [{ op: "line", x: [1, 0], y: [0, 1] }],
	}));
	assert.equal(view.layers.get("grid").commands.length, 2);
});

test("append with a mismatched z or kind is rejected", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope({ layer: { id: "grid", action: "replace", z: -10, kind: "grid" } }));
	assert.throws(() => view.apply(envelope({
		layer: { id: "grid", action: "append", z: 5 }, commands: [],
	})), EnvelopeFormatError);
	assert.throws(() => view.apply(envelope({
		layer: { id: "grid", action: "append", kind: "points" }, commands: [],
	})), EnvelopeFormatError);
});

test("replace keeps the creation slot among z-ties", () => {
	const view = new CornishView(makeSVG());
	view.apply(envelope({ layer: { id: "a", action: "replace" } }));
	view.apply(envelope({ layer: { id: "b", action: "replace" } }));
	const original = view.layers.get("a").seq;
	view.apply(envelope({ layer: { id: "a", action: "replace" }, commands: [] }));
	assert.equal(view.layers.get("a").seq, original);
});

test("delete removes the layer group from the DOM", () => {
	const svg = makeSVG();
	const view = new CornishView(svg);
	view.apply(envelope({ layer: { id: "grid", action: "replace" } }));
	assert.equal(svg.querySelectorAll('[data-layer-id="grid"]').length, 1);
	const body = envelope({ layer: { id: "grid", action: "delete" } });
	delete body.commands;
	delete body.canvas;
	view.apply(body);
	assert.equal(svg.querySelectorAll('[data-layer-id="grid"]').length, 0);
});

test("layers render sorted by z then creation order", () => {
	const svg = makeSVG();
	const view = new CornishView(svg);
	view.apply({
		cornish_view: "1.0", type: "draw_stream", canvas: CANVAS,
		layers: [
			{ id: "top", z: 5, commands: [] },
			{ id: "back", z: -20, commands: [] },
			{ id: "mid", commands: [] },
		],
	});
	const order = svg.descendants()
		.filter((node) => node.nodeName === "g" && node.getAttribute("data-layer-id"))
		.map((node) => node.getAttribute("data-layer-id"));
	assert.deepEqual(order, ["back", "mid", "top"]);
});

test("a layer with no kind omits data-kind entirely", () => {
	const svg = makeSVG();
	new CornishView(svg).apply(envelope({ layer: { id: "bare", action: "replace" } }));
	const group = svg.querySelectorAll('[data-layer-id="bare"]')[0];
	assert.equal(group.getAttribute("data-kind"), null);
});

// ---------------------------------------------------------------------------
// the frozen fallbacks (§3.3) — behaviours, with a diagnostic
// ---------------------------------------------------------------------------

const FALLBACK_COMMANDS = [
	{ op: "attr", prim: "line", colour: 999, style: 9 },
	{ op: "line", x: [0.1, 0.9], y: [0.1, 0.9] },
	{ op: "attr", prim: "text", font: 42 },
	{ op: "text", text: "x", x: 0.5, y: 0.5, just: "CC", up: [0, 1] },
	{ op: "mark", x: [0.5], y: [0.5], marker: 99 },
	{ op: "hologram", intensity: 3 },
];

test("every out-of-range value degrades visibly and warns", () => {
	const svg = makeSVG();
	const view = new CornishView(svg);
	const warnings = captureWarnings(() => {
		view.apply(envelope({ layer: { id: "fallbacks", action: "replace" }, commands: FALLBACK_COMMANDS }));
	});

	const joined = warnings.join("\n");
	for (const fragment of ["colour index 999", "line style 9", "font index 42", "marker code 99", "unknown command op"]) {
		assert.ok(joined.includes(fragment), `no console.warn mentioned ${fragment}\n${joined}`);
	}

	const nodes = svg.descendants();
	const polyline = nodes.find((node) => node.nodeName === "polyline");
	assert.ok(polyline, "the line still rendered");
	assert.equal(polyline.getAttribute("stroke"), "#000000", "colour 999 -> foreground");
	assert.equal(polyline.getAttribute("stroke-dasharray"), null, "style 9 -> solid");
	assert.ok(nodes.some((node) => node.nodeName === "text"), "the text still rendered");
	assert.ok(nodes.some((node) => node.nodeName === "circle"), "marker 99 -> a dot");
});

// ---------------------------------------------------------------------------
// the frozen mapping tables (§4.2) — must agree with the Python writer
// ---------------------------------------------------------------------------

function renderCommands(commands, canvas = CANVAS) {
	const svg = makeSVG();
	new CornishView(svg).apply(envelope({
		canvas, layer: { id: "l", action: "replace" }, commands,
	}));
	return svg;
}

test("line styles map to the frozen dasharrays", () => {
	for (const [style, dashed] of [[1, false], [2, true], [3, true], [4, true], [5, true]]) {
		const svg = renderCommands([
			{ op: "attr", prim: "line", style },
			{ op: "line", x: [0.1, 0.9], y: [0.1, 0.9] },
		]);
		const polyline = svg.descendants().find((node) => node.nodeName === "polyline");
		assert.equal(polyline.getAttribute("stroke-dasharray") !== null, dashed, `style ${style}`);
	}
});

test("font indices map to the frozen weight/style pairs", () => {
	const expected = { 1: ["normal", "normal"], 2: ["bold", "normal"], 3: ["normal", "italic"],
		4: ["bold", "italic"], 5: ["500", "normal"], 6: ["600", "normal"] };
	for (const [font, [weight, style]] of Object.entries(expected)) {
		const svg = renderCommands([
			{ op: "attr", prim: "text", font: Number(font) },
			{ op: "text", text: "x", x: 0.5, y: 0.5, just: "CC", up: [0, 1] },
		]);
		const text = svg.descendants().find((node) => node.nodeName === "text");
		assert.equal(text.getAttribute("font-weight"), weight, `font ${font} weight`);
		assert.equal(text.getAttribute("font-style"), style, `font ${font} style`);
	}
});

test("justification maps to the frozen SVG attributes", () => {
	const cases = [["TL", "start", "hanging"], ["CC", "middle", "middle"],
		["BR", "end", "alphabetic"], ["ML", "start", "text-after-edge"]];
	for (const [just, anchor, baseline] of cases) {
		const svg = renderCommands([{ op: "text", text: "x", x: 0.5, y: 0.5, just, up: [0, 1] }]);
		const text = svg.descendants().find((node) => node.nodeName === "text");
		assert.equal(text.getAttribute("text-anchor"), anchor, just);
		assert.equal(text.getAttribute("dominant-baseline"), baseline, just);
	}
});

test("every frozen marker code produces geometry", () => {
	const expected = { 1: "circle", 2: "path", 3: "path", 4: "circle", 5: "path",
		6: "circle", 7: "polygon", 8: "polygon", 9: "polygon", 10: "polygon" };
	for (const [code, tag] of Object.entries(expected)) {
		const svg = renderCommands([{ op: "mark", x: [0.5], y: [0.5], marker: Number(code) }]);
		assert.ok(svg.descendants().some((node) => node.nodeName === tag), `marker ${code} -> ${tag}`);
	}
});

test("horizontal text is unrotated and vertical text rotates -90", () => {
	let svg = renderCommands([{ op: "text", text: "x", x: 0.5, y: 0.5, just: "CC", up: [0, 1] }]);
	assert.equal(svg.descendants().find((n) => n.nodeName === "text").getAttribute("transform"), null);
	svg = renderCommands([{ op: "text", text: "x", x: 0.5, y: 0.5, just: "CC", up: [-1, 0] }]);
	const transform = svg.descendants().find((n) => n.nodeName === "text").getAttribute("transform");
	assert.ok(transform && transform.startsWith("rotate(-90"), transform);
});

test("the viewBox follows the §4.3 coordinate mapping", () => {
	const svg = renderCommands([], { gbox: [0, 0, 1, 0.5], aspect: 0.5, margin: 0.15, width: 640 });
	const [, , width, height] = svg.getAttribute("viewBox").split(" ").map(Number);
	assert.equal(width, 640);
	// (Gh + 2m) / (Gw + 2m) with Gw = 1, Gh = 0.5, m = 0.15
	assert.ok(Math.abs(height / width - 0.8 / 1.3) < 1e-6, `${height}/${width}`);
});

test("grid lines and marks clip to the gbox while text does not", () => {
	const svg = renderCommands([
		{ op: "line", x: [0.1, 0.9], y: [0.1, 0.9] },
		{ op: "text", text: "x", x: 0.5, y: -0.1, just: "CC", up: [0, 1] },
	]);
	const nodes = svg.descendants();
	assert.ok(nodes.find((n) => n.nodeName === "polyline").getAttribute("clip-path"));
	assert.equal(nodes.find((n) => n.nodeName === "text").getAttribute("clip-path"), null);
});

test("a remote image href is rejected", () => {
	assert.throws(() => renderCommands([
		{ op: "image", x1: 0, y1: 0, x2: 1, y2: 1, href: "https://example.org/x.png" },
	]), EnvelopeFormatError);
});

test("clear removes everything", () => {
	const svg = makeSVG();
	const view = new CornishView(svg);
	view.apply(envelope({ commands: [{ op: "line", x: [0, 1], y: [0, 1] }] }));
	view.clear();
	assert.equal(svg.childNodes.length, 0);
	assert.equal(view.layers.size, 0);
});
