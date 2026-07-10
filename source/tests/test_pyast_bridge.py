
'''
Unit tests for the private bridge module (SPEC-04A §9, tiers T1–T5).

The parametrized tables below compile the SPEC-04A §2 decision table directly;
test ids follow the spec's row ids (P1–P20, K1–K6).
'''

import numpy as np
import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord
import starlink.Ast as Ast

from cornish import _pyast_bridge as bridge
from cornish import ASTCircle, ASTSkyFrame, ASTFrame as CornishFrame
from cornish.mapping.frame.sky_frame import ASTICRSFrame


# ---------------------------------------------------------------------------
# fixtures: the six frame kinds (§2.1), built once per session
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def K1():
	''' sky frame '''
	return Ast.SkyFrame("System=ICRS")

@pytest.fixture(scope="session")
def K2():
	''' basic frame '''
	return Ast.Frame(2)

@pytest.fixture(scope="session")
def K3(K1, K2):
	''' frameset, current frame sky '''
	fs = Ast.FrameSet(K2)
	fs.addframe(1, Ast.UnitMap(2), K1)
	assert fs.getframe(Ast.CURRENT).isaskyframe()
	return fs

@pytest.fixture(scope="session")
def K4(K1, K2):
	''' frameset, current frame basic '''
	fs = Ast.FrameSet(K1)
	fs.addframe(1, Ast.UnitMap(2), K2)
	assert not fs.getframe(Ast.CURRENT).isaskyframe()
	return fs

@pytest.fixture(scope="session")
def K5(K1):
	''' region on a sky frame '''
	return Ast.Circle(K1, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(2.0)])

@pytest.fixture(scope="session")
def K6(K2):
	''' region on a basic frame '''
	return Ast.Circle(K2, 1, [50.0, 50.0], [10.0])


SKY_KINDS = ["K1", "K3", "K5"]
BASIC_KINDS = ["K2", "K4", "K6"]
ALL_KINDS = SKY_KINDS + BASIC_KINDS

@pytest.fixture
def frames(K1, K2, K3, K4, K5, K6):
	return {"K1": K1, "K2": K2, "K3": K3, "K4": K4, "K5": K5, "K6": K6}


# reference point and array (degrees)
REF_POINT = np.array([10.0, 20.0])
REF_ARRAY = np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]])
REF_RAD_SINGLE = np.deg2rad(REF_POINT).reshape(2, 1)      # (2,1)
REF_RAD_PARALLEL = np.deg2rad(REF_ARRAY.T)                # (2,3)


def galactic_equivalent_of_icrs(ra_deg, dec_deg):
	''' A galactic-frame SkyCoord equal to the given ICRS position. '''
	return SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs").galactic


# ---------------------------------------------------------------------------
# T1 — decision-table matrix
# ---------------------------------------------------------------------------

class TestT1DecisionTable:

	# ---- input builders; one per P row that has a defined sky-kind result ----

	def check_sky(self, points, kind_frame, expected, **kwargs):
		result = bridge.to_frame_units(points, kind_frame, **kwargs)
		assert result.dtype == np.float64
		assert result.flags['C_CONTIGUOUS']
		np.testing.assert_allclose(result, expected, rtol=0, atol=1e-12)

	@pytest.mark.parametrize("kind", SKY_KINDS)
	def test_P1_scalar_skycoord_sky(self, frames, kind):
		sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg, frame="icrs")
		self.check_sky(sc, frames[kind], REF_RAD_SINGLE)

	@pytest.mark.parametrize("kind", BASIC_KINDS)
	def test_P1_scalar_skycoord_basic(self, frames, kind):
		sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg, frame="icrs")
		with pytest.raises(ValueError):
			bridge.to_frame_units(sc, frames[kind])

	@pytest.mark.parametrize("kind", SKY_KINDS)
	def test_P1g_galactic_skycoord_honored(self, frames, kind):
		''' a galactic SkyCoord into an ICRS frame is converted, not assumed '''
		sc = galactic_equivalent_of_icrs(10.0, 20.0)
		result = bridge.to_frame_units(sc, frames[kind])
		np.testing.assert_allclose(result, REF_RAD_SINGLE, atol=1e-10)

	@pytest.mark.parametrize("kind", SKY_KINDS)
	def test_P2_array_skycoord(self, frames, kind):
		sc = SkyCoord(ra=REF_ARRAY[:, 0] * u.deg, dec=REF_ARRAY[:, 1] * u.deg, frame="icrs")
		self.check_sky(sc, frames[kind], REF_RAD_PARALLEL)

	def test_P2_multidim_skycoord_valueerror(self, K1):
		sc = SkyCoord(ra=[[10, 30], [50, 70]] * u.deg, dec=[[20, 40], [60, 80]] * u.deg)
		with pytest.raises(ValueError):
			bridge.to_frame_units(sc, K1)

	@pytest.mark.parametrize("kind", BASIC_KINDS)
	def test_P2_array_skycoord_basic(self, frames, kind):
		sc = SkyCoord(ra=REF_ARRAY[:, 0] * u.deg, dec=REF_ARRAY[:, 1] * u.deg, frame="icrs")
		with pytest.raises(ValueError):
			bridge.to_frame_units(sc, frames[kind])

	def test_P3_mixed_frame_skycoord_sequence(self, K1):
		''' each element converts independently — no shared-frame requirement '''
		sc1 = SkyCoord(ra=10 * u.deg, dec=20 * u.deg, frame="icrs")
		sc2 = galactic_equivalent_of_icrs(30.0, 40.0)
		result = bridge.to_frame_units([sc1, sc2], K1)
		expected = np.deg2rad([[10.0, 30.0], [20.0, 40.0]])
		np.testing.assert_allclose(result, expected, atol=1e-10)

	def test_P3_basic_valueerror(self, K2):
		sc1 = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
		with pytest.raises(ValueError):
			bridge.to_frame_units([sc1, sc1], K2)

	@pytest.mark.parametrize("kind", SKY_KINDS)
	def test_P4_scalar_quantity_pair(self, frames, kind):
		''' units unified across the sequence (V7) '''
		self.check_sky((10 * u.deg, 1200 * u.arcmin), frames[kind], REF_RAD_SINGLE)

	def test_P4_basic_valueerror(self, K2):
		with pytest.raises(ValueError):
			bridge.to_frame_units((10 * u.deg, 20 * u.deg), K2)

	def test_P5_parallel_quantity_arrays(self, K1):
		points = (np.r_[10., 30., 50.] * u.deg, np.r_[20., 40., 60.] * u.deg)
		result = bridge.to_frame_units(points, K1, parallel_axes=True)
		np.testing.assert_allclose(result, REF_RAD_PARALLEL, atol=1e-12)

	def test_P5_basic_valueerror(self, K2):
		points = (np.r_[10., 30., 50.] * u.deg, np.r_[20., 40., 60.] * u.deg)
		with pytest.raises(ValueError):
			bridge.to_frame_units(points, K2, parallel_axes=True)

	@pytest.mark.parametrize("kind", SKY_KINDS)
	def test_P6_rank1_quantity(self, frames, kind):
		self.check_sky([10, 20] * u.deg, frames[kind], REF_RAD_SINGLE)

	@pytest.mark.parametrize("kind", BASIC_KINDS)
	def test_P6_basic_valueerror(self, frames, kind):
		with pytest.raises(ValueError):
			bridge.to_frame_units([10, 20] * u.deg, frames[kind])

	def test_P7_pairs_quantity(self, K1):
		self.check_sky(REF_ARRAY * u.deg, K1, REF_RAD_PARALLEL)

	def test_P8_parallel_quantity(self, K1):
		self.check_sky(REF_ARRAY.T * u.deg, K1, REF_RAD_PARALLEL)

	def test_P9_square_quantity_reads_as_pairs(self, K1):
		''' D3: (2,2) is pairs — column 0 is the first point '''
		q = np.array([[10.0, 20.0], [30.0, 40.0]]) * u.deg
		result = bridge.to_frame_units(q, K1)
		np.testing.assert_allclose(result[:, 0], np.deg2rad([10.0, 20.0]), atol=1e-12)
		np.testing.assert_allclose(result[:, 1], np.deg2rad([30.0, 40.0]), atol=1e-12)

	@pytest.mark.parametrize("kind", ALL_KINDS)
	def test_P10_bare_single_point(self, frames, kind):
		expected = REF_RAD_SINGLE if kind in SKY_KINDS else REF_POINT.reshape(2, 1)
		result = bridge.to_frame_units([10.0, 20.0], frames[kind])
		np.testing.assert_allclose(result, expected, atol=1e-12)

	@pytest.mark.parametrize("kind", ALL_KINDS)
	def test_P11_bare_pairs(self, frames, kind):
		expected = REF_RAD_PARALLEL if kind in SKY_KINDS else REF_ARRAY.T
		result = bridge.to_frame_units(REF_ARRAY, frames[kind])
		np.testing.assert_allclose(result, expected, atol=1e-12)

	@pytest.mark.parametrize("kind", ALL_KINDS)
	def test_P12_bare_parallel(self, frames, kind):
		expected = REF_RAD_PARALLEL if kind in SKY_KINDS else REF_ARRAY.T
		result = bridge.to_frame_units(REF_ARRAY.T, frames[kind])
		np.testing.assert_allclose(result, expected, atol=1e-12)

	@pytest.mark.parametrize("sky", [True, False])
	def test_P13_square_bare(self, K1, K2, sky):
		frame = K1 if sky else K2
		square = np.array([[10.0, 20.0], [30.0, 40.0]])
		def maybe_rad(x):
			return np.deg2rad(x) if sky else np.asarray(x, dtype=float)
		# no flag: pairs
		result = bridge.to_frame_units(square, frame)
		np.testing.assert_allclose(result[:, 0], maybe_rad([10.0, 20.0]), atol=1e-12)
		# parallel_axes=True: parallel
		result = bridge.to_frame_units(square, frame, parallel_axes=True)
		np.testing.assert_allclose(result[:, 0], maybe_rad([10.0, 30.0]), atol=1e-12)
		# parallel_axes=False: pairs
		result = bridge.to_frame_units(square, frame, parallel_axes=False)
		np.testing.assert_allclose(result[:, 0], maybe_rad([10.0, 20.0]), atol=1e-12)

	@pytest.mark.parametrize("kind", ALL_KINDS)
	def test_P14_string(self, frames, kind):
		s = "((10,20),(30,40),(50,60))"
		expected = REF_RAD_PARALLEL if kind in SKY_KINDS else REF_ARRAY.T
		result = bridge.to_frame_units(s, frames[kind])
		np.testing.assert_allclose(result, expected, atol=1e-12)

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	def test_P15_nonangular_quantity(self, frames, kind):
		with pytest.raises(ValueError):
			bridge.to_frame_units([10, 20] * u.m, frames[kind])

	def test_P15_chained_from_unitconversionerror(self, K1):
		with pytest.raises(ValueError) as excinfo:
			bridge.to_frame_units([10, 20] * u.m, K1)
		assert isinstance(excinfo.value.__cause__, u.UnitConversionError)

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	def test_P16_dimensionless_quantity(self, frames, kind):
		with pytest.raises(ValueError):
			bridge.to_frame_units(u.Quantity([10., 20.]), frames[kind])

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	@pytest.mark.parametrize("empty", [[], np.empty((0, 2)), np.empty((2, 0))],
	                         ids=["list", "0x2", "2x0"])
	def test_P17_empty(self, frames, kind, empty):
		with pytest.raises(ValueError):
			bridge.to_frame_units(empty, frames[kind])

	def test_P17_empty_skycoord_and_quantity(self, K1):
		with pytest.raises(ValueError):
			bridge.to_frame_units(SkyCoord(ra=[] * u.deg, dec=[] * u.deg), K1)
		with pytest.raises(ValueError):
			bridge.to_frame_units(u.Quantity([]), K1)

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	@pytest.mark.parametrize("bad_value", [[np.nan, 20], [10, np.inf], [Ast.BAD, 20]],
	                         ids=["nan", "inf", "astbad"])
	def test_P18_nonfinite(self, frames, kind, bad_value):
		with pytest.raises(ValueError):
			bridge.to_frame_units(bad_value, frames[kind])

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	@pytest.mark.parametrize("bad_shape", [np.ones((3, 4)), np.float64(5), np.ones((2, 2, 2))],
	                         ids=["3x4", "rank0", "rank3"])
	def test_P19_bad_shapes(self, frames, kind, bad_shape):
		with pytest.raises(ValueError):
			bridge.to_frame_units(bad_shape, frames[kind])

	@pytest.mark.parametrize("kind", ["K1", "K2"])
	@pytest.mark.parametrize("bad_input", [{"ra": 1}, None, [True, False], ["a", "b"]],
	                         ids=["dict", "none", "bools", "strs"])
	def test_P20_typeerror(self, frames, kind, bad_input):
		with pytest.raises(TypeError):
			bridge.to_frame_units(bad_input, frames[kind])

	def test_frame_not_framelike(self):
		with pytest.raises(TypeError):
			bridge.to_frame_units([10.0, 20.0], Ast.UnitMap(2))

	# ---- frame-kind welds ----

	def test_weld_K3_equals_K1(self, K1, K3):
		np.testing.assert_array_equal(bridge.to_frame_units([10.0, 20.0], K3),
		                              bridge.to_frame_units([10.0, 20.0], K1))

	def test_weld_K4_equals_K2(self, K2, K4):
		np.testing.assert_array_equal(bridge.to_frame_units([10.0, 20.0], K4),
		                              bridge.to_frame_units([10.0, 20.0], K2))

	def test_weld_K5_equals_K1(self, K1, K5):
		np.testing.assert_array_equal(bridge.to_frame_units([10.0, 20.0], K5),
		                              bridge.to_frame_units([10.0, 20.0], K1))

	def test_weld_K6_equals_K2(self, K2, K6):
		np.testing.assert_array_equal(bridge.to_frame_units([10.0, 20.0], K6),
		                              bridge.to_frame_units([10.0, 20.0], K2))

	def test_footgun_pin(self, K3, K4):
		'''
		is_sky must dispatch on the CURRENT frame; also pin the raw pyast footgun
		(FrameSet.isaskyframe() is False even when the current frame is sky) so a
		future pyast behavior change is NOTICED rather than silently double-guarded.
		'''
		assert bridge.is_sky(K3) is True
		assert bridge.is_sky(K4) is False
		assert K3.isaskyframe() is False

	def test_cornish_wrapper_unwrapping(self, K1):
		''' cornish-wrapper variants pin _unwrap on a subset (P1, P10) '''
		wrapper = ASTSkyFrame(ast_object=K1)
		sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
		np.testing.assert_allclose(bridge.to_frame_units(sc, wrapper), REF_RAD_SINGLE, atol=1e-12)
		np.testing.assert_allclose(bridge.to_frame_units([10.0, 20.0], wrapper),
		                           REF_RAD_SINGLE, atol=1e-12)
		basic_wrapper = CornishFrame(naxes=2)
		np.testing.assert_allclose(bridge.to_frame_units([10.0, 20.0], basic_wrapper),
		                           REF_POINT.reshape(2, 1), atol=1e-12)

	def test_naxes3_frame(self):
		''' naxes > 2 frames get the same shape rules with naxes substituted (V3) '''
		f3 = Ast.Frame(3)
		result = bridge.to_frame_units([1.0, 2.0, 3.0], f3)
		np.testing.assert_array_equal(result, [[1.0], [2.0], [3.0]])
		result = bridge.to_frame_units(np.arange(12.0).reshape(4, 3), f3)
		assert result.shape == (3, 4)
		with pytest.raises(ValueError):
			bridge.to_frame_units([1.0, 2.0], f3)


# ---------------------------------------------------------------------------
# T2 — squeeze / parallel_axes / single-point-form edges + dispatch-order pins
# ---------------------------------------------------------------------------

class TestT2Edges:

	def test_squeeze_P1(self, K1):
		sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
		result = bridge.to_frame_units(sc, K1, squeeze=True)
		assert result.shape == (2,)
		np.testing.assert_allclose(result, np.deg2rad(REF_POINT), atol=1e-12)

	def test_squeeze_P10(self, K1):
		result = bridge.to_frame_units([10.0, 20.0], K1, squeeze=True)
		assert result.shape == (2,)

	def test_squeeze_multipoint_valueerror(self, K1):
		with pytest.raises(ValueError):
			bridge.to_frame_units(REF_ARRAY, K1, squeeze=True)

	@pytest.mark.parametrize("points", [
		pytest.param(SkyCoord(ra=10 * u.deg, dec=20 * u.deg), id="P1"),
		pytest.param((10 * u.deg, 20 * u.arcmin), id="P4"),
		pytest.param("((10,20),(30,40))", id="P14"),
		pytest.param([10.0, 20.0], id="P10"),
	])
	def test_parallel_axes_invalid_forms(self, K1, points):
		with pytest.raises(ValueError):
			bridge.to_frame_units(points, K1, parallel_axes=True)

	def test_parallel_axes_contradictions(self, K1):
		with pytest.raises(ValueError):
			bridge.to_frame_units(np.ones((3, 2)), K1, parallel_axes=True)
		with pytest.raises(ValueError):
			bridge.to_frame_units(np.ones((2, 3)), K1, parallel_axes=False)

	def test_is_single_point_true_forms(self, K1):
		assert bridge.is_single_point(SkyCoord(ra=10 * u.deg, dec=20 * u.deg), K1) is True
		assert bridge.is_single_point([10., 20.], K1) is True
		assert bridge.is_single_point([10, 20] * u.deg, K1) is True
		assert bridge.is_single_point((10 * u.deg, 20 * u.arcmin), K1) is True
		assert bridge.is_single_point("(10, 20)", K1) is True

	def test_is_single_point_false_forms(self, K1):
		sc_arr = SkyCoord(ra=[10, 30] * u.deg, dec=[20, 40] * u.deg)
		assert bridge.is_single_point(sc_arr, K1) is False
		assert bridge.is_single_point([[10., 20.]], K1) is False          # (1,2) array form
		assert bridge.is_single_point(np.ones((2, 1)), K1) is False
		assert bridge.is_single_point(REF_ARRAY, K1) is False
		assert bridge.is_single_point("((10,20),(30,40))", K1) is False
		sc = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
		assert bridge.is_single_point([sc], K1) is False                  # sequence form

	def test_is_single_point_invalid_input(self, K1):
		with pytest.raises(TypeError):
			bridge.is_single_point({"ra": 1}, K1)

	def test_dispatch_order_bare_list_not_swallowed_by_P16(self, K1):
		'''
		Round-1 regression class: u.Quantity([10.0, 20.0]) SUCCEEDS on a bare list
		(dimensionless), so coercion-success dispatch would swallow every bare
		array into the P16 ValueError. Bare degrees must convert.
		'''
		result = bridge.to_frame_units([10.0, 20.0], K1)
		np.testing.assert_allclose(result, REF_RAD_SINGLE, atol=1e-12)

	def test_dispatch_order_string_reaches_parser(self, K1):
		''' string input must reach parse_points_string, not die in Quantity coercion '''
		result = bridge.to_frame_units("(10, 20)", K1)
		np.testing.assert_allclose(result, REF_RAD_SINGLE, atol=1e-12)

	def test_incompatible_quantity_sequence_valueerror(self, K1):
		''' step-5 stacking failure semantics (V17): ValueError chained from UnitConversionError '''
		with pytest.raises(ValueError) as excinfo:
			bridge.to_frame_units((10 * u.deg, 5 * u.m), K1)
		assert isinstance(excinfo.value.__cause__, u.UnitConversionError)


# ---------------------------------------------------------------------------
# T3 — distance pair
# ---------------------------------------------------------------------------

class TestT3DistancePair:

	def test_float_on_sky_is_degrees(self, K1):
		assert bridge.to_frame_distance(1.0, K1) == pytest.approx(np.deg2rad(1.0))

	def test_quantity_on_sky(self, K1):
		assert bridge.to_frame_distance(1 * u.arcmin, K1) == pytest.approx(np.deg2rad(1 / 60))

	def test_float_on_basic_native(self, K2):
		assert bridge.to_frame_distance(1.0, K2) == 1.0

	def test_int_accepted(self, K1):
		assert bridge.to_frame_distance(2, K1) == pytest.approx(np.deg2rad(2.0))

	def test_quantity_on_basic_valueerror(self, K2):
		with pytest.raises(ValueError):
			bridge.to_frame_distance(1 * u.deg, K2)

	def test_nonangular_quantity_valueerror(self, K1):
		with pytest.raises(ValueError):
			bridge.to_frame_distance(1 * u.m, K1)

	@pytest.mark.parametrize("value", [np.nan, np.inf, Ast.BAD])
	def test_nonfinite_valueerror(self, K1, value):
		with pytest.raises(ValueError):
			bridge.to_frame_distance(value, K1)

	def test_bool_typeerror(self, K1):
		with pytest.raises(TypeError):
			bridge.to_frame_distance(True, K1)

	def test_nonscalar_quantity_valueerror(self, K1):
		with pytest.raises(ValueError):
			bridge.to_frame_distance([1, 2] * u.deg, K1)

	def test_from_frame_distance_sky(self, K1):
		assert bridge.from_frame_distance(np.deg2rad(1.0), K1) == pytest.approx(1.0)

	def test_from_frame_distance_basic(self, K2):
		assert bridge.from_frame_distance(1.5, K2) == 1.5

	@pytest.mark.parametrize("value", [np.nan, np.inf, Ast.BAD])
	def test_from_frame_distance_nonfinite(self, K1, value):
		with pytest.raises(ValueError):
			bridge.from_frame_distance(value, K1)


# ---------------------------------------------------------------------------
# T4 — is_sky / current_frame_of (§3 table + liveness pins)
# ---------------------------------------------------------------------------

class TestT4EffectiveFrame:

	def test_region_current_frame(self, K5):
		frame = bridge.current_frame_of(K5)
		assert frame.isaskyframe()
		assert bridge.is_sky(K5) is True
		assert K5.isaskyframe() is False  # the V1 footgun, pinned

	def test_region_on_basic(self, K6):
		assert bridge.is_sky(K6) is False
		assert not bridge.current_frame_of(K6).isaskyframe()

	def test_plain_frame_identity(self, K1):
		assert bridge.current_frame_of(K1) is K1

	def test_nonframe_typeerror(self):
		with pytest.raises(TypeError):
			bridge.current_frame_of(Ast.UnitMap(2))
		with pytest.raises(TypeError):
			bridge.is_sky(Ast.UnitMap(2))
		with pytest.raises(TypeError):
			bridge.current_frame_of(None)

	def test_cornish_wrapper(self):
		circle = ASTCircle(center=[30, 45], radius=2.0)
		assert bridge.is_sky(circle) is True

	def test_frameset_live_pointer(self, K2):
		'''
		Documented live behavior (V2): mutating the copy=False return of a
		FrameSet DOES propagate. This test exists so a pyast change to
		copy-semantics is caught.
		'''
		fs = Ast.FrameSet(K2)
		fs.addframe(1, Ast.UnitMap(2), Ast.SkyFrame("System=ICRS"))
		live = bridge.current_frame_of(fs)
		assert fs.getframe(Ast.CURRENT).get("System") == "ICRS"
		live.set("System=Galactic")
		assert fs.getframe(Ast.CURRENT).get("System") == "GALACTIC"

	def test_region_frame_is_copy(self, K1):
		region = Ast.Circle(K1, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(2.0)])
		frame = bridge.current_frame_of(region)
		frame.set("System=Galactic")
		assert region.getregionframe().get("System") == "ICRS"

	def test_copy_true_never_propagates(self, K2):
		fs = Ast.FrameSet(K2)
		fs.addframe(1, Ast.UnitMap(2), Ast.SkyFrame("System=ICRS"))
		independent = bridge.current_frame_of(fs, copy=True)
		independent.set("System=Galactic")
		assert fs.getframe(Ast.CURRENT).get("System") == "ICRS"


# ---------------------------------------------------------------------------
# T5 — dump_value
# ---------------------------------------------------------------------------

class TestT5DumpValue:

	@pytest.fixture(scope="class")
	def sky(self):
		return Ast.SkyFrame("System=ICRS")

	@pytest.fixture(scope="class")
	def circle_pair(self, sky):
		c1 = Ast.Circle(sky, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(1.0)])
		c2 = Ast.Circle(sky, 1, np.deg2rad([32.0, 44.0]), [np.deg2rad(1.0)])
		return c1, c2

	def test_operator_labels_and_weld(self, circle_pair):
		''' the dump label is 'Operator' (1=AND, 2=OR) — welded to the Ast constants '''
		c1, c2 = circle_pair
		cmp_or = Ast.CmpRegion(c1, c2, Ast.OR)
		cmp_and = Ast.CmpRegion(c1, c2, Ast.AND)
		assert bridge.dump_value(cmp_or, "Operator") == "2"
		assert bridge.dump_value(cmp_and, "Operator") == "1"
		assert int(bridge.dump_value(cmp_or, "Operator")) == Ast.OR
		assert int(bridge.dump_value(cmp_and, "Operator")) == Ast.AND

	def test_unc_block_readback(self, sky):
		'''
		The Unc block is Channel-readable; AST internally reshapes the stored
		uncertainty region, so the read-back size is only approximate (~6% at
		this dec-20 centre — the SPEC-04A §7.3 verified basis; never compare
		exactly).
		'''
		from cornish.channel.channel_io import ListSource
		unc = Ast.Circle(sky, 1, np.deg2rad([10.0, 20.0]), [np.deg2rad(0.001)])
		c = Ast.Circle(sky, 1, np.deg2rad([10.0, 20.0]), [np.deg2rad(1.0)], unc)  # POSITIONAL
		block = bridge.dump_value(c, "Unc")
		read_back = Ast.Channel(ListSource(block), None).read()
		radius_deg = np.rad2deg(read_back.circlepars()[1])
		assert 0.0009 <= radius_deg <= 0.0012

	def test_unc_missing_without_positional_unc(self, sky):
		'''
		Unc appears only when an uncertainty region was positionally supplied:
		a raw pyast region with no unc argument runs on the AST-internal default
		and KeyError is the honest answer.
		'''
		c = Ast.Circle(sky, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(1.0)])
		with pytest.raises(KeyError):
			bridge.dump_value(c, "Unc")

	def test_include_defaults(self, sky):
		''' System is a '#' default line on a region dump '''
		c = Ast.Circle(sky, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(1.0)])
		assert bridge.dump_value(c, "System", include_defaults=True) == "ICRS"
		with pytest.raises(KeyError):
			bridge.dump_value(c, "System", include_defaults=False)

	def test_missing_component_keyerror(self, sky):
		c = Ast.Circle(sky, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(1.0)])
		with pytest.raises(KeyError):
			bridge.dump_value(c, "NoSuchComponent")

	def test_nested_component_not_visible(self, sky):
		''' depth rule: components of nested records (e.g. the SkyFrame inside a region) are not depth-1 '''
		c = Ast.Circle(sky, 1, np.deg2rad([30.0, 45.0]), [np.deg2rad(1.0)])
		with pytest.raises(KeyError):
			bridge.dump_value(c, "SkyTol")

	def test_typeerror_for_non_ast(self):
		with pytest.raises(TypeError):
			bridge.dump_value("not an ast object", "System")

	def test_quoted_scalar_unquoted(self, sky):
		''' surrounding double quotes are removed from scalar values '''
		assert bridge.dump_value(sky, "System") == "ICRS"
