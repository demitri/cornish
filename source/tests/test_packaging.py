"""Packaging check (cornish issue #3): the built wheel must contain every
file under the cornish package tree (subpackages and any data files)."""
import pathlib
import subprocess
import sys
import zipfile

SOURCE_DIR = pathlib.Path(__file__).resolve().parents[1]


def test_wheel_contains_whole_package_tree(tmp_path):
	subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
					"-w", str(tmp_path), str(SOURCE_DIR)], check=True)
	(wheel,) = tmp_path.glob("cornish-*.whl")
	shipped = set(zipfile.ZipFile(wheel).namelist())
	tree = {p.relative_to(SOURCE_DIR).as_posix()
			for p in (SOURCE_DIR / "cornish").rglob("*")
			if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
	assert tree - shipped == set()
