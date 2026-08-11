from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUILD_DIR = ROOT / "build"
if BUILD_DIR.exists():
    shutil.rmtree(BUILD_DIR)

def run(cmd: list[str], **kwargs) -> None:
    print("+", " ".join(cmd))
    subprocess.check_call(cmd, **kwargs)

def findCmake()-> str:
    cmake = shutil.which('cmake')
    if cmake:
        return cmake
    raise FileNotFoundError('cmake not found on PATH.\nIn your conda env run:\n conda install -c conda-forge cmake ninja\nThen open a *new* terminal, conda activate science, and retry.')

def findNinja()-> str | None:
    return shutil.which('ninja')

def ensurePyBind11()-> str:
    try:
        import pybind11
    except ImportError:
        print('Installing pybind11 for this Python...')
        run([sys.executable, '-m', 'pip', 'install', "pybind11"])
        import pybind11
    return pybind11.get_cmake_dir()

def findBuiltModule(buildDirectory)-> Path | None:
    preferred = list((buildDirectory / 'python').glob('engine.*'))
    if preferred:
        return preferred[0]
    candidates = list(buildDirectory.rglob('engine*.pyd')) + list(buildDirectory.rglob('engine*.so'))
    candidates = [c for c in candidates if c.suffix.lower() in {'.pyd', '.so'}]
    return candidates[0] if candidates else None


def build(clean: bool = False)-> int:
    print(f'Using Python {sys.version}')
    print(f'Executable: {sys.executable}')
    if "3.10" not in sys.version:
        print("WARNING: this is not Python 3.10. \nThe .pyd will only load in the same major.minor version.")
    cmake = findCmake()
    ninja = findNinja()
    print(f"cmake: {cmake}")
    print(f"ninja: {ninja}" if ninja else "ninja not on PATH — CMake will pick the default generator.")
    pybind11Directory = ensurePyBind11()
    print(f"pybind11 cmake: {pybind11Directory}")
    if clean and BUILD_DIR.exists():
        print(f"Cleaning {BUILD_DIR}")
        shutil.rmtree(BUILD_DIR)

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    cmakeConfiguration = ['cmake', '-S', str(ROOT), '-B', str(BUILD_DIR), f'-DPython_EXECUTABLE={sys.executable}', f'-Dpybind11_DIR={pybind11Directory}', '-DCMAKE_BUILD_TYPE=Release']
    cmakeConfiguration = cmakeConfiguration + ["-G", "Ninja"] if ninja else cmakeConfiguration
    run(cmakeConfiguration)
    run(['cmake', '--build', str(BUILD_DIR), '--config', 'Release'])
    module = findBuiltModule(BUILD_DIR)
    if module is None:
        print("ERROR: built engine module not found under", BUILD_DIR)
        print("Check the CMake/build log for compile or link errors.")
        return 1

    destination = ROOT / module.name
    shutil.copy2(module, destination)
    print(f'Installed: {destination}')
    sys.path.insert(0, str(ROOT))
    if 'engine' in sys.modules:
        del sys.modules['engine']
    import engine
    print('engine.__file__ =', engine.__file__)
    print('has SimulatorParameters:', hasattr(engine, 'SimulatorParameters'))
    if not hasattr(engine, 'SimulatorParameters'):
        print('ERROR: module loaded but SimulatorParameters is missing (bindings incomplete?)')
        return 1
    print('OK — engine is ready for this Python.')
    return 0

if __name__ == "__main__":
    raise SystemExit(build())