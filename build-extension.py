import os
import shlex
import shutil
import subprocess
import zipfile
from pathlib import Path

""" 
Note: this build script is invoked by Poetry during the build process
source: https://python-poetry.org/docs/building-extension-modules/#maturin-build-script
"""

def maturin(*args):
    subprocess.call(["maturin", *list(args)])


def build():
    build_dir = Path(__file__).parent.joinpath("build")
    build_dir.mkdir(parents=True, exist_ok=True)

    wheels_dir = Path(__file__).parent.joinpath("target/wheels")
    if wheels_dir.exists():
        shutil.rmtree(wheels_dir)

    cargo_args = []
    if os.getenv("MATURIN_BUILD_ARGS"):
        cargo_args = shlex.split(os.getenv("MATURIN_BUILD_ARGS", ""))

    maturin("build", "-r", *cargo_args)

    # We won't use the wheel built by maturin directly since
    # we want Poetry to build it, but we need to retrieve the
    # compiled extensions from the maturin wheel.
    wheel = next(iter(wheels_dir.glob("*.whl")))
    with zipfile.ZipFile(wheel.as_posix()) as whl:
        whl.extractall(wheels_dir.as_posix())

        for extension in wheels_dir.rglob("**/*.so"):
            shutil.copyfile(extension, build_dir.joinpath(extension.name))

    # shutil.rmtree(wheels_dir)


if __name__ == "__main__":
    build()