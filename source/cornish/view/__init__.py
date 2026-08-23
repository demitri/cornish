
'''
``cornish.view`` — serving plots to a browser (SPEC-09A §1, §5).

Import direction is binding: **``cornish.view`` imports ``cornish.plot``, never
the reverse** at module level. That separation is what keeps ``render_svg`` a
pure call usable in any web framework (Flask pattern 1, the recommended
production shape) with no transport dependency dragged along.

Within this package, only :mod:`~cornish.view.server` and
:mod:`~cornish.view.blueprint` touch transport, and ``blueprint`` is the only
module that imports Flask — behind a guarded import with an actionable message,
so a missing Flask is never a mysterious traceback.

``flask_blueprint`` is re-exported lazily here for the same reason: importing
``cornish.view`` must not require Flask.
'''

from ..plot.svg import render_svg
from .registry import PlotRegistry, Subscription

__all__ = ['render_svg', 'PlotRegistry', 'Subscription', 'flask_blueprint']


def flask_blueprint(*args, **kwargs):
	'''
	Lazy passthrough to :func:`cornish.view.blueprint.flask_blueprint`, so that
	importing this package does not require Flask.
	'''
	from .blueprint import flask_blueprint as build
	return build(*args, **kwargs)
