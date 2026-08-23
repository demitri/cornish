
'''
``PlotRegistry`` — the in-process store of plot documents (SPEC-09A §1.1, §3.4).

Transport-free by design: the stdlib daemon and the Flask blueprint share this
verbatim, which is what makes SPEC-09's "the blueprint shares 100% of its
handler code with the standalone daemon" true rather than aspirational.

Two rules shape everything here:

* **A plot exists if and only if it has received a canvas.** §3.1 requires a
  canvas on the envelope that creates a plot, so no code path can ever observe
  a canvas-less document.
* **Nothing is silently ignored.** An append or delete addressing a layer that
  is not there is an error, not a no-op: an append is a statement about
  existing content, and auto-creating would mask a lost earlier envelope.
'''

import queue
import threading
from typing import Optional

from ..exc import EnvelopeError
from ..plot.envelope import (ID_PATTERN, SCHEMA_VERSION, LayerAction, validate_document,
                             validate_envelope)
from ..plot.svg import render_svg

__all__ = ['PlotRegistry', 'Subscription', 'SUBSCRIBER_QUEUE_LIMIT']

#: How many events a lagging subscriber may fall behind before its connection is
#: dropped. Dropping forces a clean snapshot reconnect, which is strictly better
#: than silently discarding events and leaving the viewer subtly wrong forever.
SUBSCRIBER_QUEUE_LIMIT = 256


class Subscription:
	'''
	A queue-backed subscription to one plot's event stream.

	Events are ``(event_name, payload)`` pairs, where *event_name* is one of
	``snapshot``, ``envelope`` or ``plot_deleted``.
	'''

	def __init__(self, plot_id: str) -> None:
		self.plot_id = plot_id
		self.queue = queue.Queue(maxsize=SUBSCRIBER_QUEUE_LIMIT)
		#: set when the queue overflowed; the server drops the connection
		self.overflowed = False
		#: True while this subscriber still owes an initial snapshot. It stays
		#: True for a subscriber that arrived before the plot existed — the
		#: DEFERRED snapshot case (§5.3).
		self.snapshot_pending = True

	def _offer(self, event: str, payload: dict) -> None:
		if self.overflowed:
			return                      # already doomed; the server is about to drop it
		try:
			self.queue.put_nowait((event, payload))
		except queue.Full:
			self.overflowed = True


class PlotRegistry:
	'''
	Thread-safe store of draw-stream documents, keyed by plot id.
	'''

	def __init__(self) -> None:
		self._lock = threading.RLock()
		self._plots = {}
		self._subscribers = {}

	# ------------------------------------------------------------------
	# reads
	# ------------------------------------------------------------------

	def plot_ids(self) -> list:
		''' The ids of plots that exist — i.e. that have received a canvas. '''
		with self._lock:
			return sorted(self._plots)

	def document(self, plot_id: str) -> dict:
		'''
		The full draw-stream document, layers in render order (§3.5).

		:raises KeyError: if no such plot exists
		'''
		with self._lock:
			return self._document_locked(plot_id)

	def revision(self, plot_id: str) -> int:
		'''
		The plot's current revision.

		:raises KeyError: if no such plot exists
		'''
		with self._lock:
			return self._require_plot(plot_id)["revision"]

	def svg(self, plot_id: str, *, width: int = 640) -> str:
		'''
		The plot rendered to SVG.

		:raises KeyError: if no such plot exists
		'''
		return render_svg(self.document(plot_id), width=width)

	def _require_plot(self, plot_id: str) -> dict:
		plot = self._plots.get(plot_id)
		if plot is None:
			raise KeyError(f"No plot '{plot_id}'.")
		return plot

	def _document_locked(self, plot_id: str) -> dict:
		plot = self._require_plot(plot_id)
		ordered = sorted(plot["layers"].items(), key=lambda item: (item[1]["z"], item[1]["seq"]))
		document = {
			"cornish_view": SCHEMA_VERSION,
			"type": "draw_stream",
			"plot_id": plot_id,
			"canvas": plot["canvas"],
			"layers": [{"id": layer_id, "kind": layer["kind"], "z": layer["z"],
			            "commands": layer["commands"]}
			           if layer["kind"] is not None else
			           {"id": layer_id, "z": layer["z"], "commands": layer["commands"]}
			           for layer_id, layer in ordered],
		}
		return document

	# ------------------------------------------------------------------
	# writes
	# ------------------------------------------------------------------

	def apply(self, plot_id: str, envelope_or_document: dict) -> int:
		'''
		Apply one envelope, or one full document, to *plot_id*.

		A **document** atomically replaces the entire plot (canvas, then every
		layer in order) — this is what makes ``curl -d @saved.json`` and
		``show()``'s bootstrap work, since a ``.json`` export is a document.

		Shape dispatch is by the §3.1 marker key: a body carrying ``layers`` is
		a document, anything else is an envelope (which covers the canvas-only
		form). Each validator then rejects the other shape's marker key, so a
		hybrid body fails whichever way it is dispatched.

		:returns: the plot's new revision
		:raises EnvelopeError: if the body is not a valid envelope/document
		:raises ValueError: if the body carries a plot_id that is not *plot_id*
		:raises KeyError: for append/delete on an unknown layer, or any envelope for an unknown plot
		'''
		if not ID_PATTERN.match(str(plot_id)):
			raise ValueError(f"'{plot_id}' is not a legal plot id (must match {ID_PATTERN.pattern}).")
		if not isinstance(envelope_or_document, dict):
			raise EnvelopeError(f"/: expected a JSON object; got {type(envelope_or_document).__name__}.")

		is_document = "layers" in envelope_or_document
		# validation precedes the plot_id check, matching the frozen §5.2
		# failure-mode precedence (parse/validation -> mismatch -> state)
		if is_document:
			validate_document(envelope_or_document)
		else:
			validate_envelope(envelope_or_document)

		body_plot_id = envelope_or_document.get("plot_id")
		if body_plot_id is not None and body_plot_id != plot_id:
			raise ValueError(f"The body's plot_id '{body_plot_id}' does not match the "
			                 f"addressed plot '{plot_id}'.")

		with self._lock:
			if is_document:
				revision, event, payload = self._apply_document(plot_id, envelope_or_document)
			else:
				revision, event, payload = self._apply_envelope(plot_id, envelope_or_document)
			self._fan_out(plot_id, event, payload)
			return revision

	def _apply_document(self, plot_id: str, document: dict):
		''' Atomic whole-plot replace. '''
		previous = self._plots.get(plot_id)
		revision = (previous["revision"] + 1) if previous else 1

		layers = {}
		for index, layer in enumerate(document["layers"], start=1):
			layers[layer["id"]] = {
				"kind": layer.get("kind"),
				"z": int(layer.get("z", 0)),
				"commands": list(layer["commands"]),
				"seq": index,
			}
		self._plots[plot_id] = {
			"canvas": dict(document["canvas"]),
			"layers": layers,
			"sequence": len(layers),
			"revision": revision,
		}
		# deltas cannot express "everything you knew is gone", so a document
		# apply fans out as a fresh snapshot (§5.3)
		return revision, "snapshot", self._document_locked(plot_id)

	def _apply_envelope(self, plot_id: str, envelope: dict):
		canvas = envelope.get("canvas")
		plot = self._plots.get(plot_id)

		if plot is None:
			if canvas is None:
				# a plot exists iff it has received a canvas, so an envelope
				# without one cannot create it
				raise KeyError(f"No plot '{plot_id}': the envelope that creates a plot must "
				               "carry a canvas.")
			plot = {"canvas": dict(canvas), "layers": {}, "sequence": 0, "revision": 0}
			self._plots[plot_id] = plot
		elif canvas is not None:
			# a later canvas fully replaces the stored one and dirties every layer
			plot["canvas"] = dict(canvas)

		layer_spec = envelope.get("layer")
		if layer_spec is not None:
			self._apply_layer(plot, layer_spec, envelope.get("commands"))

		plot["revision"] += 1
		return plot["revision"], "envelope", envelope

	def _apply_layer(self, plot: dict, layer_spec: dict, commands) -> None:
		layer_id = layer_spec["id"]
		action = LayerAction(layer_spec["action"])
		existing = plot["layers"].get(layer_id)

		if action == LayerAction.DELETE:
			if existing is None:
				raise KeyError(f"Cannot delete layer '{layer_id}': it does not exist.")
			del plot["layers"][layer_id]
			return

		if action == LayerAction.APPEND:
			if existing is None:
				# auto-creating would silently mask a lost earlier envelope
				raise KeyError(f"Cannot append to layer '{layer_id}': it does not exist. "
				               "An append is a statement about existing content.")
			# on append, z and kind must be absent or agree — a layer changing
			# type or depth mid-life is a producer bug worth surfacing (§3.4)
			if "z" in layer_spec and int(layer_spec["z"]) != existing["z"]:
				raise EnvelopeError(f"/layer/z: {layer_spec['z']} does not match the stored z "
				                    f"({existing['z']}) for layer '{layer_id}'; on append it must "
				                    "be absent or equal.")
			if "kind" in layer_spec and layer_spec["kind"] != existing["kind"]:
				raise EnvelopeError(f"/layer/kind: '{layer_spec['kind']}' does not match the "
				                    f"stored kind ('{existing['kind']}') for layer '{layer_id}'; "
				                    "on append it must be absent or equal.")
			existing["commands"].extend(commands)
			return

		# replace: create, or overwrite content/kind/z while RETAINING the
		# creation sequence, so a restyled layer keeps its slot among z-ties
		if existing is None:
			plot["sequence"] += 1
			sequence = plot["sequence"]
		else:
			sequence = existing["seq"]
		plot["layers"][layer_id] = {
			"kind": layer_spec.get("kind"),
			"z": int(layer_spec.get("z", 0)),
			"commands": list(commands),
			"seq": sequence,
		}

	def delete_plot(self, plot_id: str) -> None:
		'''
		Remove a plot and tell its viewers.

		:raises KeyError: if no such plot exists
		'''
		with self._lock:
			self._require_plot(plot_id)
			del self._plots[plot_id]
			self._fan_out(plot_id, "plot_deleted", {})

	# ------------------------------------------------------------------
	# subscriptions
	# ------------------------------------------------------------------

	def subscribe(self, plot_id: str) -> Subscription:
		'''
		Subscribe to a plot's events.

		Subscribing to a **not-yet-POSTed** id creates a subscription *slot*,
		not a plot: the id does not appear in :meth:`plot_ids`,
		:meth:`document` still raises ``KeyError``, and the snapshot is
		DEFERRED until the first canvas-bearing envelope creates the plot. A
		viewer opened before the first POST is a supported workflow.
		'''
		subscription = Subscription(plot_id)
		with self._lock:
			self._subscribers.setdefault(plot_id, []).append(subscription)
			if plot_id in self._plots:
				subscription._offer("snapshot", self._document_locked(plot_id))
				subscription.snapshot_pending = False
		return subscription

	def unsubscribe(self, subscription: Subscription) -> None:
		''' Drop a subscription. Unsubscribing twice is harmless. '''
		with self._lock:
			listeners = self._subscribers.get(subscription.plot_id, [])
			if subscription in listeners:
				listeners.remove(subscription)
			if not listeners:
				self._subscribers.pop(subscription.plot_id, None)

	def subscriber_count(self, plot_id: str) -> int:
		with self._lock:
			return len(self._subscribers.get(plot_id, []))

	def _fan_out(self, plot_id: str, event: str, payload: dict) -> None:
		'''
		Deliver an event to every subscriber, **after** the write has been
		committed — a slow viewer can lag, never corrupt.

		One frozen exception (§5.3): a subscriber still awaiting its deferred
		snapshot receives the plot-creating envelope as a ``snapshot`` built
		from the just-committed document, not as an ``envelope``. The creating
		envelope IS the whole document by §3.1, so the content is equivalent,
		and the viewer's snapshot-then-deltas contract then holds for every
		subscriber regardless of when it subscribed.
		'''
		listeners = self._subscribers.get(plot_id)
		if not listeners:
			return
		document = None
		for subscription in list(listeners):
			if event == "plot_deleted":
				subscription._offer(event, payload)
				subscription.snapshot_pending = True     # a later plot needs a fresh snapshot
				continue
			if subscription.snapshot_pending:
				if document is None:
					document = self._document_locked(plot_id)
				subscription._offer("snapshot", document)
				subscription.snapshot_pending = False
			else:
				subscription._offer(event, payload)
