"""Packaging check (cornish issue #3): the built wheel must contain every
file under the cornish package tree, and every module in the unpacked wheel
must import."""
import pathlib
import subprocess
import sys
import zipfile

SOURCE_DIR = pathlib.Path(__file__).resolve().parents[1]

# optional third-party packages a module may legitimately be missing
OPTIONAL_THIRD_PARTY = {"cornish.plot.matplotlib": "matplotlib"}

IMPORT_ALL = '''
import importlib, pkgutil, sys
sys.path.insert(0, sys.argv[1])
import cornish
failures = []
for m in pkgutil.walk_packages(cornish.__path__, "cornish."):
	try:
		importlib.import_module(m.name)
	except ModuleNotFoundError as e:
		if e.name == {optional!r}.get(m.name):
			continue
		failures.append((m.name, repr(e)))
	except Exception as e:
		failures.append((m.name, repr(e)))
assert cornish.__file__.startswith(sys.argv[1]), cornish.__file__
if failures:
	sys.exit("import failures: %r" % failures)
'''


def test_wheel_contains_and_imports_whole_package(tmp_path):
	subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
					"-w", str(tmp_path), str(SOURCE_DIR)], check=True)
	(wheel,) = tmp_path.glob("cornish-*.whl")
	zf = zipfile.ZipFile(wheel)
	tree = {p.relative_to(SOURCE_DIR).as_posix()
			for p in (SOURCE_DIR / "cornish").rglob("*")
			if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
	assert tree - set(zf.namelist()) == set()

	site = tmp_path / "site"
	zf.extractall(site)
	subprocess.run([sys.executable, "-I", "-c", IMPORT_ALL.format(optional=OPTIONAL_THIRD_PARTY), str(site)],
				   check=True, cwd=tmp_path)
