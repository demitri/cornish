/**
 * cornish-view.js — the browser replay renderer (SPEC-09A §4).
 *
 * Dependency-free ESM, no build step, served as one static file.
 *
 * This is a REPLAY client. AST already did all the layout server-side: the
 * renderer does no text measurement, no label placement and no projection
 * maths. It implements the same §4.2–4.4 mapping tables as the Python SVG
 * writer — one table, two implementations, cross-checked by the shared
 * conformance fixtures (§7.6).
 */

const SVG_NS = "http://www.w3.org/2000/svg";
const SCHEMA_MAJOR = 1;
const CLIP_PATH_ID = "cornish-gbox";

/** Line style → stroke-dasharray, in multiples of the stroke width (§4.2). */
const DASH_PATTERNS = {
	1: [],
	2: [8, 4],
	3: [8, 4, 2, 4],
	4: [2, 4],
	5: [8, 4, 2, 4, 2, 4],
};

/** Font index → [font-weight, font-style] (§3.3). */
const FONT_STYLES = {
	1: ["normal", "normal"],
	2: ["bold", "normal"],
	3: ["normal", "italic"],
	4: ["bold", "italic"],
	5: ["500", "normal"],
	6: ["600", "normal"],
};

const BASELINES = { T: "hanging", C: "middle", B: "alphabetic", M: "text-after-edge" };
const ANCHORS = { L: "start", C: "middle", R: "end" };

const MARKER_RADIUS_UNIT = 0.006;
const STROKE_WIDTH_UNIT = 0.002;
const MINIMUM_STROKE_PX = 0.5;

const ATTR_DEFAULTS = { colour: 1, width: 1.0, style: 1, size: 1.0, font: 1 };

const DEFAULT_PALETTE = [
	"rgba(0,0,0,0)", "#000000", "#d62728", "#2ca02c",
	"#1f77b4", "#17becf", "#9467bd", "#bcbd22",
	"#ff7f0e", "#8c564b", "#e377c2", "#7f7f7f",
	"#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5",
];

/** Raised for any malformed envelope or document. */
export class EnvelopeFormatError extends Error {
	constructor(message) {
		super(message);
		this.name = "EnvelopeFormatError";
	}
}

/**
 * RFC 8785 (JCS) canonical serialization (§3.6).
 *
 * The reference JS implementation: recursive sorted-key stringify. ES number
 * formatting is native here — which is exactly why the Python side implements
 * the ECMAScript Number::toString algorithm rather than using its own repr.
 * Keys sort by UTF-16 code units, which is what Array#sort does for strings.
 */
export function canonicalJSON(value) {
	if (value === null) return "null";
	const type = typeof value;
	if (type === "boolean") return value ? "true" : "false";
	if (type === "number") {
		if (!Number.isFinite(value)) {
			throw new EnvelopeFormatError(`Non-finite number ${value} cannot appear in an envelope.`);
		}
		return String(value === 0 ? 0 : value); // collapses -0 to 0, per JCS
	}
	if (type === "string") return JSON.stringify(value);
	if (Array.isArray(value)) return "[" + value.map(canonicalJSON).join(",") + "]";
	if (type === "object") {
		const keys = Object.keys(value).sort();
		return "{" + keys.map((key) => JSON.stringify(key) + ":" + canonicalJSON(value[key])).join(",") + "}";
	}
	throw new EnvelopeFormatError(`Values of type ${type} cannot appear in an envelope.`);
}

function formatNumber(value) {
	let text = Number(value).toFixed(3);
	if (text.indexOf(".") !== -1) {
		text = text.replace(/0+$/, "").replace(/\.$/, "");
	}
	return text === "" || text === "-0" ? "0" : text;
}

function diagnose(message) {
	// §3.7: a fallback is announced, never silent. The JS channel is console.warn;
	// render_svg's is an SVG comment plus a logged warning. One definition, two
	// channels, so independent consumers cannot invent incompatible ones.
	if (typeof console !== "undefined" && console.warn) {
		console.warn("cornish-view: " + message);
	}
}

export class CornishView {
	/**
	 * @param {SVGSVGElement} svgElement the target; CornishView owns its content
	 * @param {object} [options]
	 * @param {string[]} [options.palette] base-palette theme override (literal CSS
	 *        colours carried in a stream are never remapped)
	 * @param {string} [options.background] override canvas.background
	 */
	constructor(svgElement, options = {}) {
		if (!svgElement) throw new EnvelopeFormatError("CornishView needs a target <svg> element.");
		this.svg = svgElement;
		this.options = options;
		this.document = svgElement.ownerDocument || (typeof document !== "undefined" ? document : null);
		if (!this.document) throw new EnvelopeFormatError("No DOM document is available.");
		this.clear();
	}

	/** Remove all layers and canvas state. */
	clear() {
		this.canvas = null;
		this.layers = new Map(); // id -> {kind, z, commands, seq}
		this.sequence = 0;
		this.layerElements = new Map();
		while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);
	}

	/**
	 * Apply one envelope OR a full document.
	 *
	 * Shape dispatch replicates §3.1's exactly-one-of rule INDEPENDENTLY of the
	 * Python validators (which this cannot call): a body carrying both `layer`
	 * and `layers` is rejected either way round.
	 */
	apply(body) {
		if (body === null || typeof body !== "object" || Array.isArray(body)) {
			throw new EnvelopeFormatError("An envelope must be a JSON object.");
		}
		this._checkVersion(body);
		if (body.type !== "draw_stream") {
			throw new EnvelopeFormatError(`Expected type 'draw_stream'; got '${body.type}'.`);
		}

		const hasLayer = Object.prototype.hasOwnProperty.call(body, "layer");
		const hasLayers = Object.prototype.hasOwnProperty.call(body, "layers");
		if (hasLayer && hasLayers) {
			throw new EnvelopeFormatError(
				"A body carries exactly one of 'layer' (envelope) or 'layers' (document), never both.");
		}

		if (hasLayers) {
			this._applyDocument(body);
		} else {
			this._applyEnvelope(body, hasLayer);
		}
		return this;
	}

	_checkVersion(body) {
		const version = body.cornish_view;
		if (typeof version !== "string" || !/^\d+\.\d+$/.test(version)) {
			throw new EnvelopeFormatError(`'${version}' is not a '<major>.<minor>' version string.`);
		}
		const major = parseInt(version.split(".")[0], 10);
		if (major !== SCHEMA_MAJOR) {
			// unknown MINOR is accepted (§3.7); a MAJOR changes existing semantics
			throw new EnvelopeFormatError(`stream requires cornish-view ${major}.x; this renderer implements ${SCHEMA_MAJOR}.x.`);
		}
	}

	_applyDocument(body) {
		if (!body.canvas) throw new EnvelopeFormatError("A document requires a 'canvas'.");
		if (!Array.isArray(body.layers)) throw new EnvelopeFormatError("'layers' must be an array.");
		this.clear();
		this.canvas = body.canvas;
		const seen = new Set();
		body.layers.forEach((layer, index) => {
			if (!layer || typeof layer.id !== "string") {
				throw new EnvelopeFormatError(`layers[${index}] needs an 'id'.`);
			}
			if (seen.has(layer.id)) {
				throw new EnvelopeFormatError(`Duplicate layer id '${layer.id}'.`);
			}
			if (!Array.isArray(layer.commands)) {
				throw new EnvelopeFormatError(`layers[${index}] needs a 'commands' array.`);
			}
			seen.add(layer.id);
			this.layers.set(layer.id, {
				kind: layer.kind,
				z: layer.z === undefined ? 0 : layer.z,
				commands: layer.commands.slice(),
				seq: ++this.sequence,
			});
		});
		this._renderAll();
	}

	_applyEnvelope(body, hasLayer) {
		const canvasChanged = Boolean(body.canvas);
		if (canvasChanged) this.canvas = body.canvas;

		if (!hasLayer) {
			// the canvas-only form: plot creation / re-theming
			if (!canvasChanged) {
				throw new EnvelopeFormatError(
					"An envelope needs a 'layer', unless it is the canvas-only form.");
			}
			if (Object.prototype.hasOwnProperty.call(body, "commands")) {
				throw new EnvelopeFormatError("A canvas-only envelope must not carry 'commands'.");
			}
			this._renderAll();
			return;
		}

		const spec = body.layer;
		if (!spec || typeof spec.id !== "string") {
			throw new EnvelopeFormatError("'layer' needs an 'id'.");
		}
		const action = spec.action;
		if (["replace", "append", "delete"].indexOf(action) === -1) {
			throw new EnvelopeFormatError(`'${action}' is not a legal layer action.`);
		}

		const existing = this.layers.get(spec.id);
		if (action === "delete") {
			if (!existing) throw new EnvelopeFormatError(`Cannot delete layer '${spec.id}': it does not exist.`);
			this.layers.delete(spec.id);
			if (canvasChanged) { this._renderAll(); return; }
			this._removeLayerElement(spec.id);
			return;
		}

		if (!Array.isArray(body.commands)) {
			throw new EnvelopeFormatError(`A '${action}' envelope needs a 'commands' array.`);
		}

		if (action === "append") {
			if (!existing) {
				throw new EnvelopeFormatError(
					`Cannot append to layer '${spec.id}': it does not exist. An append is a statement about existing content.`);
			}
			if (spec.z !== undefined && spec.z !== existing.z) {
				throw new EnvelopeFormatError(`layer.z ${spec.z} does not match the stored z (${existing.z}) on append.`);
			}
			if (spec.kind !== undefined && spec.kind !== existing.kind) {
				throw new EnvelopeFormatError(`layer.kind '${spec.kind}' does not match the stored kind on append.`);
			}
			existing.commands = existing.commands.concat(body.commands);
		} else {
			this.layers.set(spec.id, {
				kind: spec.kind,
				z: spec.z === undefined ? 0 : spec.z,
				commands: body.commands.slice(),
				// a restyled layer keeps its slot among z-ties
				seq: existing ? existing.seq : ++this.sequence,
			});
		}

		// a canvas change dirties every layer; otherwise only this one is rebuilt
		if (canvasChanged) this._renderAll();
		else this._renderLayerIncremental(spec.id);
	}

	/** The serialized current SVG (export / save-as). */
	toSVGString() {
		return serializeElement(this.svg);
	}

	// -- geometry (§4.3) ------------------------------------------------

	_geometry() {
		const canvas = this.canvas || {};
		const gbox = canvas.gbox || [0, 0, 1, 1];
		const [x1, y1, x2, y2] = gbox.map(Number);
		const gw = x2 - x1;
		const gh = y2 - y1;
		const margin = (canvas.margin === undefined ? 0.15 : Number(canvas.margin)) * Math.max(gw, gh);
		const width = canvas.width === undefined ? 640 : Number(canvas.width);
		const scale = width / (gw + 2 * margin);
		const height = scale * (gh + 2 * margin);
		return { x1, y1, x2, y2, margin, scale, width, height };
	}

	_palette() {
		if (this.options.palette) return this.options.palette;
		if (this.canvas && this.canvas.palette) return this.canvas.palette;
		return DEFAULT_PALETTE;
	}

	_colour(value) {
		if (typeof value === "string") return value;
		const palette = this._palette();
		const index = Number(value);
		if (index >= 0 && index < palette.length) return palette[index];
		diagnose(`colour index ${index} is outside the ${palette.length}-entry palette; rendering as foreground (index 1)`);
		return palette.length > 1 ? palette[1] : "#000000";
	}

	// -- rendering ------------------------------------------------------

	_renderAll() {
		while (this.svg.firstChild) this.svg.removeChild(this.svg.firstChild);
		this.layerElements = new Map();
		if (!this.canvas) return;

		const geometry = this._geometry();
		this.svg.setAttribute("viewBox", `0 0 ${formatNumber(geometry.width)} ${formatNumber(geometry.height)}`);
		this.svg.setAttribute("width", "100%");
		this.svg.setAttribute("preserveAspectRatio", "xMidYMid meet");

		const defs = this._element("defs");
		const clip = this._element("clipPath");
		clip.setAttribute("id", CLIP_PATH_ID);
		const rect = this._element("rect");
		rect.setAttribute("x", formatNumber(this._px(geometry.x1, geometry)));
		rect.setAttribute("y", formatNumber(this._py(geometry.y2, geometry)));
		rect.setAttribute("width", formatNumber(geometry.scale * (geometry.x2 - geometry.x1)));
		rect.setAttribute("height", formatNumber(geometry.scale * (geometry.y2 - geometry.y1)));
		clip.appendChild(rect);
		defs.appendChild(clip);
		this.svg.appendChild(defs);

		const background = this.options.background !== undefined
			? this.options.background
			: (this.canvas.background === undefined ? null : this.canvas.background);
		if (background !== null && background !== undefined) {
			const fill = this._element("rect");
			fill.setAttribute("x", "0");
			fill.setAttribute("y", "0");
			fill.setAttribute("width", formatNumber(geometry.width));
			fill.setAttribute("height", formatNumber(geometry.height));
			fill.setAttribute("fill", String(background));
			this.svg.appendChild(fill);
		}

		for (const id of this._orderedLayerIds()) {
			const group = this._buildLayer(id, geometry);
			this.svg.appendChild(group);
			this.layerElements.set(id, group);
		}
	}

	_orderedLayerIds() {
		return Array.from(this.layers.entries())
			.sort((a, b) => (a[1].z - b[1].z) || (a[1].seq - b[1].seq))
			.map((entry) => entry[0]);
	}

	_removeLayerElement(id) {
		const element = this.layerElements.get(id);
		if (element && element.parentNode) element.parentNode.removeChild(element);
		this.layerElements.delete(id);
	}

	_renderLayerIncremental(id) {
		if (!this.canvas) { this._renderAll(); return; }
		const geometry = this._geometry();
		const rebuilt = this._buildLayer(id, geometry);
		const previous = this.layerElements.get(id);
		if (previous && previous.parentNode) {
			previous.parentNode.removeChild(previous);
		}
		// insert in render order rather than always appending, so z is honoured
		const order = this._orderedLayerIds();
		const position = order.indexOf(id);
		let inserted = false;
		for (let index = position + 1; index < order.length; index += 1) {
			const following = this.layerElements.get(order[index]);
			if (following && following.parentNode === this.svg) {
				this.svg.insertBefore(rebuilt, following);
				inserted = true;
				break;
			}
		}
		if (!inserted) this.svg.appendChild(rebuilt);
		this.layerElements.set(id, rebuilt);
	}

	_element(name) {
		return this.document.createElementNS(SVG_NS, name);
	}

	_px(x, geometry) { return geometry.scale * (Number(x) - geometry.x1 + geometry.margin); }
	_py(y, geometry) { return geometry.height - geometry.scale * (Number(y) - geometry.y1 + geometry.margin); }

	_buildLayer(id, geometry) {
		const layer = this.layers.get(id);
		const group = this._element("g");
		// §3.5: NO id attribute — the legal id charset contains '.' and '-',
		// which break '#...' selectors. Lookup is by attribute selector.
		group.setAttribute("data-layer-id", id);
		if (layer.kind !== undefined && layer.kind !== null) {
			group.setAttribute("data-kind", String(layer.kind));
		}

		// attr state is per-layer and resets at every layer head (§3.4)
		const state = {
			line: Object.assign({}, ATTR_DEFAULTS),
			mark: Object.assign({}, ATTR_DEFAULTS),
			text: Object.assign({}, ATTR_DEFAULTS),
		};

		for (const command of layer.commands) {
			this._renderCommand(command, state, group, geometry);
		}
		return group;
	}

	_strokePx(widthAttr, geometry) {
		return Math.max(STROKE_WIDTH_UNIT * Number(widthAttr) * geometry.scale, MINIMUM_STROKE_PX);
	}

	_renderCommand(command, state, group, geometry) {
		const op = command.op;
		if (op === "attr") {
			const target = state[command.prim];
			if (!target) { diagnose(`unknown attr prim '${command.prim}' ignored`); return; }
			for (const key of ["colour", "width", "style", "size", "font"]) {
				if (command[key] !== undefined) target[key] = command[key];
			}
			return;
		}
		if (op === "line") { group.appendChild(this._buildLine(command, state.line, geometry)); return; }
		if (op === "mark") { this._buildMarks(command, state.mark, group, geometry); return; }
		if (op === "text") { group.appendChild(this._buildText(command, state.text, geometry)); return; }
		if (op === "image") { group.appendChild(this._buildImage(command, geometry)); return; }
		// §3.7: skipped WITH a diagnostic, never silently and never fatally
		diagnose(`unknown command op '${op}' skipped (this consumer implements schema 1.0)`);
	}

	_dash(style, strokePx) {
		let pattern = DASH_PATTERNS[Number(style)];
		if (pattern === undefined) {
			diagnose(`line style ${style} is outside the frozen 1-5 table; rendering solid`);
			pattern = DASH_PATTERNS[1];
		}
		if (!pattern.length) return null;
		return pattern.map((value) => formatNumber(value * strokePx)).join(",");
	}

	_buildLine(command, state, geometry) {
		const element = this._element("polyline");
		const points = [];
		for (let index = 0; index < command.x.length; index += 1) {
			points.push(`${formatNumber(this._px(command.x[index], geometry))},${formatNumber(this._py(command.y[index], geometry))}`);
		}
		const strokePx = this._strokePx(state.width, geometry);
		element.setAttribute("points", points.join(" "));
		element.setAttribute("fill", "none");
		element.setAttribute("stroke", this._colour(state.colour));
		element.setAttribute("stroke-width", formatNumber(strokePx));
		const dash = this._dash(state.style, strokePx);
		if (dash) element.setAttribute("stroke-dasharray", dash);
		element.setAttribute("clip-path", `url(#${CLIP_PATH_ID})`);
		return element;
	}

	_buildMarks(command, state, group, geometry) {
		const colour = this._colour(state.colour);
		const strokePx = this._strokePx(state.width, geometry);
		const radius = MARKER_RADIUS_UNIT * Number(state.size) * geometry.scale;
		let code = Number(command.marker);
		if (!(code >= 1 && code <= 10)) {
			diagnose(`marker code ${code} is outside the frozen 1-10 table; rendering as a dot`);
			code = 6;
		}
		const wrapper = this._element("g");
		wrapper.setAttribute("clip-path", `url(#${CLIP_PATH_ID})`);
		for (let index = 0; index < command.x.length; index += 1) {
			const cx = this._px(command.x[index], geometry);
			const cy = this._py(command.y[index], geometry);
			for (const element of this._markerElements(code, cx, cy, radius, colour, strokePx)) {
				wrapper.appendChild(element);
			}
		}
		group.appendChild(wrapper);
	}

	_markerElements(code, cx, cy, radius, colour, strokePx) {
		if (code === 6) {
			const dot = this._element("circle");
			dot.setAttribute("cx", formatNumber(cx));
			dot.setAttribute("cy", formatNumber(cy));
			dot.setAttribute("r", formatNumber(0.25 * radius));
			dot.setAttribute("fill", colour);
			return [dot];
		}
		if (code === 1 || code === 4) {
			const circle = this._element("circle");
			circle.setAttribute("cx", formatNumber(cx));
			circle.setAttribute("cy", formatNumber(cy));
			circle.setAttribute("r", formatNumber((code === 1 ? 0.6 : 1.0) * radius));
			circle.setAttribute("fill", "none");
			circle.setAttribute("stroke", colour);
			circle.setAttribute("stroke-width", formatNumber(strokePx));
			return [circle];
		}
		const segments = { 2: [0, 90], 5: [45, 135], 3: [0, 60, 120] }[code];
		if (segments) {
			const parts = [];
			for (const degrees of segments) {
				const angle = (degrees * Math.PI) / 180;
				const dx = radius * Math.cos(angle);
				const dy = radius * Math.sin(angle);
				parts.push(`M${formatNumber(cx - dx)},${formatNumber(cy - dy)}L${formatNumber(cx + dx)},${formatNumber(cy + dy)}`);
			}
			const path = this._element("path");
			path.setAttribute("d", parts.join(""));
			path.setAttribute("fill", "none");
			path.setAttribute("stroke", colour);
			path.setAttribute("stroke-width", formatNumber(strokePx));
			return [path];
		}
		// triangles 7-10, in SVG's Y-DOWN frame ("apex up" is -y)
		const apex = { 7: -90, 8: 90, 9: 180, 10: 0 }[code];
		const points = [];
		for (const offset of [0, 120, 240]) {
			const angle = ((apex + offset) * Math.PI) / 180;
			points.push(`${formatNumber(cx + radius * Math.cos(angle))},${formatNumber(cy + radius * Math.sin(angle))}`);
		}
		const polygon = this._element("polygon");
		polygon.setAttribute("points", points.join(" "));
		polygon.setAttribute("fill", "none");
		polygon.setAttribute("stroke", colour);
		polygon.setAttribute("stroke-width", formatNumber(strokePx));
		return [polygon];
	}

	_buildText(command, state, geometry) {
		const element = this._element("text");
		const charHeight = this.canvas && this.canvas.char_height !== undefined ? Number(this.canvas.char_height) : 0.02;
		const fontSize = charHeight * Number(state.size) * geometry.scale;
		const x = this._px(command.x, geometry);
		const y = this._py(command.y, geometry);
		let styles = FONT_STYLES[Number(state.font)];
		if (styles === undefined) {
			diagnose(`font index ${state.font} is outside the frozen 1-6 table; rendering as normal`);
			styles = FONT_STYLES[1];
		}
		const family = (this.canvas && this.canvas.font_family) || "DejaVu Sans";

		element.setAttribute("x", formatNumber(x));
		element.setAttribute("y", formatNumber(y));
		element.setAttribute("fill", this._colour(state.colour));
		element.setAttribute("font-family", `${family}, sans-serif`);
		element.setAttribute("font-size", formatNumber(fontSize));
		element.setAttribute("font-weight", styles[0]);
		element.setAttribute("font-style", styles[1]);
		element.setAttribute("text-anchor", ANCHORS[command.just[1]] || "middle");
		element.setAttribute("dominant-baseline", BASELINES[command.just[0]] || "middle");

		// §4.2: in SVG's Y-down frame the rotation is atan2(upx, upy)
		const angle = (Math.atan2(Number(command.up[0]), Number(command.up[1])) * 180) / Math.PI;
		if (Math.abs(angle) > 1e-9) {
			element.setAttribute("transform", `rotate(${formatNumber(angle)},${formatNumber(x)},${formatNumber(y)})`);
		}
		element.textContent = String(command.text);
		return element;
	}

	_buildImage(command, geometry) {
		if (typeof command.href !== "string" || command.href.slice(0, 5) !== "data:") {
			// envelopes are self-contained; a remote fetch would make replay
			// impure and turn a plot into a tracking vector
			throw new EnvelopeFormatError("An image href must be a 'data:' URI.");
		}
		const left = this._px(Math.min(command.x1, command.x2), geometry);
		const right = this._px(Math.max(command.x1, command.x2), geometry);
		const top = this._py(Math.max(command.y1, command.y2), geometry);
		const bottom = this._py(Math.min(command.y1, command.y2), geometry);
		const element = this._element("image");
		element.setAttribute("x", formatNumber(left));
		element.setAttribute("y", formatNumber(top));
		element.setAttribute("width", formatNumber(right - left));
		element.setAttribute("height", formatNumber(bottom - top));
		element.setAttribute("preserveAspectRatio", "none");
		element.setAttribute("clip-path", `url(#${CLIP_PATH_ID})`);
		element.setAttribute("href", command.href);
		return element;
	}
}

/** Minimal, dependency-free SVG serialization (used by toSVGString). */
function serializeElement(element) {
	if (element.nodeType === 3) return escapeText(element.nodeValue || "");
	const name = element.nodeName || element.tagName;
	const attributes = [];
	const source = element.attributes || {};
	const names = typeof source.length === "number"
		? Array.from({ length: source.length }, (_, index) => source[index].name)
		: Object.keys(source);
	for (const attributeName of names) {
		const value = element.getAttribute(attributeName);
		if (value !== null && value !== undefined) {
			attributes.push(` ${attributeName}="${escapeAttribute(String(value))}"`);
		}
	}
	const children = element.childNodes || [];
	if (!children.length) {
		const text = element.textContent;
		if (text) return `<${name}${attributes.join("")}>${escapeText(text)}</${name}>`;
		return `<${name}${attributes.join("")}/>`;
	}
	const inner = Array.from(children).map(serializeElement).join("");
	return `<${name}${attributes.join("")}>${inner}</${name}>`;
}

function escapeText(text) {
	return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function escapeAttribute(text) {
	return escapeText(text).replace(/"/g, "&quot;");
}
