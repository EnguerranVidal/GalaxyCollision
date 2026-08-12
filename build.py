from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUILD_DIR = ROOT / 'build'
if BUILD_DIR.exists():
    shutil.rmtree(BUILD_DIR)

def run(cmd: list[str], **kwargs) -> None:
    print('+', ' '.join(cmd))
    subprocess.check_call(cmd, **kwargs)

def findCmake()-> str:
    cmake = shutil.which('cmake')
    if cmake:
        return cmake
    root = Path(sys.executable).resolve().parent
    candidates = [root / 'Library' / 'bin' / 'cmake.exe',
                  root / 'Scripts' / 'cmake.exe',
                  root / 'bin' / 'cmake.exe',
                  root.parent / 'Library' / 'bin' / 'cmake.exe',
                  root.parent / 'Scripts' / 'cmake.exe']
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    raise FileNotFoundError(
        'cmake not found on PATH.\n'
        'In your conda env run:\n'
        '  conda install -c conda-forge cmake ninja\n'
        'Then open a *new* terminal, conda activate science, and retry.'
    )

def findNinja()-> str | None:
    ninja = shutil.which('ninja')
    if ninja:
        return ninja
    root = Path(sys.executable).resolve().parent
    candidates = [root / 'Library' / 'bin' / 'ninja.exe',
                  root / 'Scripts' / 'ninja.exe',
                  root.parent / 'Library' / 'bin' / 'ninja.exe']
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return None

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
    if '3.10' not in sys.version:
        print('WARNING: this is not Python 3.10. \nThe .pyd will only load in the same major.minor version.')
    cmake = findCmake()
    ninja = findNinja()
    print(f'cmake: {cmake}')
    print(f'ninja: {ninja}' if ninja else 'ninja not on PATH — CMake will pick the default generator.')
    pybind11Directory = ensurePyBind11()
    print(f'pybind11 cmake: {pybind11Directory}')
    if clean and BUILD_DIR.exists():
        print(f'Cleaning {BUILD_DIR}')
        shutil.rmtree(BUILD_DIR)

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    cmakeConfiguration = [cmake, '-S', str(ROOT), '-B', str(BUILD_DIR), f'-DPython_EXECUTABLE={sys.executable}', f'-Dpybind11_DIR={pybind11Directory}', '-DCMAKE_BUILD_TYPE=Release']
    cmakeConfiguration = cmakeConfiguration + ['-G', 'Ninja', f'-DCMAKE_MAKE_PROGRAM={ninja}'] if ninja else cmakeConfiguration
    run(cmakeConfiguration)
    buildConfiguration = [cmake, '--build', str(BUILD_DIR), '--config', 'Release']
    buildConfiguration = buildConfiguration + ['--config', 'Release'] if not ninja else buildConfiguration
    run(buildConfiguration)
    module = findBuiltModule(BUILD_DIR)
    if module is None:
        print('ERROR: built engine module not found under', BUILD_DIR)
        print('Check the CMake/build log for compile or link errors.')
        return 1

    destination = ROOT / module.name
    shutil.copy2(module, destination)
    print(f'Installed: {destination}')
    sys.path.insert(0, str(ROOT))
    if 'engine' in sys.modules:
        del sys.modules['engine']
    try:
        import engine
    except Exception as exc:
        print('ERROR: could not import engine after install:', exc)
        return 1
    print('engine.__file__ =', getattr(engine, '__file__', None))
    print('has SimulatorParameters:', hasattr(engine, 'SimulatorParameters'))
    attributes = [a for a in dir(engine) if not a.startswith("_")]
    print('exported:', ', '.join(attributes) if attributes else '(none)')
    if getattr(engine, "__file__", None) is None:
        print('ERROR: imported a namespace package (likely src/engine), not the .pyd.\nKeep the .pyd next to main.py and ensure project root is first on sys.path.')
        return 1
    if not hasattr(engine, 'SimulatorParameters'):
        print('ERROR: .pyd loaded but SimulatorParameters is missing.\nFix bindParameters.cpp / base-class bindings and rebuild.')
        return 1
    print('OK — engine is ready for this Python / GUI.')
    return 0

def main() -> int:
    parser = argparse.ArgumentParser(description='Build GalaxyCollision pybind11 engine module')
    parser.add_argument('--clean', action='store_true', help='Delete the build/ directory before configuring')
    args = parser.parse_args()
    return build(clean=args.clean)


if __name__ == '__main__':
    raise SystemExit(main())