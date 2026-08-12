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

def findVcVars() -> Path:
    candidates = [Path(r'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'),
                  Path(r'C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat'),
                  Path(r'C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat'),
                  Path(r'C:\Program Files\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat')]
    vsWhere = Path(r'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe')
    if vsWhere.is_file():
        out = subprocess.check_output([str(vsWhere), '-latest', '-products', '*', '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64', '-property', 'installationPath'], text=True).strip()
        if out:
            candidates.insert(0, Path(out) / 'VC' / 'Auxiliary' / 'Build' / 'vcvars64.bat')
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError('vcvars64.bat not found. Install VS with \"Desktop development with C++\".')

def runWithMsVc(commandList: list[str]) -> None:
    vcVars = findVcVars()
    parts = []
    for command in commandList:
        commandString = str(command)
        parts.append(f'"{commandString}"' if any(ch in commandString for ch in ' &()[]{}^=;!\'+,`~') else commandString)
    inner = " ".join(parts)
    full = f'call "{vcVars}" && {inner}'
    print("+", full)
    subprocess.check_call(full, shell=True)

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


def build(clean: bool = False, withCuda: bool = False)-> int:
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
        shutil.rmtree(BUILD_DIR)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    print('Loading MSVC via vcvars64, then configuring...')
    cmakeConfiguration = [cmake, '-S', str(ROOT), '-B', str(BUILD_DIR), f'-DPython_EXECUTABLE={sys.executable}', f'-Dpybind11_DIR={pybind11Directory}', '-DCMAKE_BUILD_TYPE=Release', f"-DGALAXY_ENABLE_CUDA={'ON' if withCuda else 'OFF'}"]
    cmakeConfiguration = cmakeConfiguration + ['-G', 'Ninja', f'-DCMAKE_MAKE_PROGRAM={ninja}'] if ninja else cmakeConfiguration
    runWithMsVc(cmakeConfiguration)
    buildConfiguration = [cmake, '--build', str(BUILD_DIR), '--config', 'Release']
    buildConfiguration = buildConfiguration + ['--config', 'Release'] if not ninja else buildConfiguration
    runWithMsVc(buildConfiguration)
    module = findBuiltModule(BUILD_DIR)
    if module is None:
        print('ERROR: no engine.pyd/.so found under', BUILD_DIR)
        return 1

    destination = ROOT / module.name
    shutil.copy2(module, destination)
    print(f'Installed: {destination}')
    sys.path.insert(0, str(ROOT))
    sys.modules.pop("engine", None)
    import engine
    print('engine.__file__ =', engine.__file__)
    ok = hasattr(engine, 'SimulatorParameters')
    print('has SimulatorParameters:', ok)
    if not ok or engine.__file__ is None:
        return 1
    print('OK — run:  python main.py')
    return 0

def main() -> int:
    parser = argparse.ArgumentParser(description='Build GalaxyCollision pybind11 engine module')
    parser.add_argument('--clean', action='store_true', help='Delete the build/ directory before configuring')
    parser.add_argument('--cuda', action='store_true', help='Enable CUDA (needs CUDA Toolkit)')
    args = parser.parse_args()
    return build(clean=args.clean, withCuda=args.cuda)


if __name__ == '__main__':
    raise SystemExit(main())