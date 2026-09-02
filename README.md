# **N-BODY SOLVER**
[![GitHub watchers](https://badgen.net/github/watchers/EnguerranVidal/GalaxyCollision/)](https://GitHub.com/EnguerranVidal/GalaxyCollision/watchers/) [![GitHub stars](https://badgen.net/github/stars/EnguerranVidal/GalaxyCollision)](https://GitHub.com/EnguerranVidal/GalaxyCollision/stargazers/)
![GitHub license](https://img.shields.io/github/license/EnguerranVidal/GalaxyCollision)

[![GitHub branches](https://badgen.net/github/branches/EnguerranVidal/GalaxyCollision)](https://github.com/EnguerranVidal/GalaxyCollision/)
[![GitHub commits](https://badgen.net/github/commits/EnguerranVidal/GalaxyCollision)](https://github.com/EnguerranVidal/GalaxyCollision/) 
![GitHub last commit](https://img.shields.io/github/last-commit/EnguerranVidal/GalaxyCollision)
![Github issus open](https://img.shields.io/github/issues-raw/EnguerranVidal/GalaxyCollision)
![Github issus closed](https://img.shields.io/github/issues-closed-raw/EnguerranVidal/GalaxyCollision)


![GitHub repo size](https://img.shields.io/github/repo-size/EnguerranVidal/GalaxyCollision)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31019/)
## SUMMARY

<div style="text-align: justify"> This repository proposes a N-Body gravity simulator with a C++/CUDA engine linked to a python PyQt5 GUI and OpenGL visualization. It supports both a Newton and a <a href="https://jheer.github.io/barnes-hut/" target="_blank">Barnes-Hut</a> force calculator with CPU and GPU backends.
The GUI and the physics engine are binded using <a href="https://pybind11.readthedocs.io/en/stable/index.html" target="_blank">pybind11</a>.  </div>

<div style="text-align: justify"> Icons by <a href="https://icons8.com/">Icons8</a>. </div>

## INSTALLATION

### Pre-Requisites

* Download [Visual Studio](https://visualstudio.microsoft.com/downloads/) (Community is enough) and select the workload: Desktop development with C++
* Install a recent [NVIDIA driver](https://www.nvidia.com/Download/index.aspx) and [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (if you want to use the GPU backend)
* Confirm the compiler is visible
```bash
nvcc --version
```

### Building the Python Environment

#### Option A — Conda (recommended on Windows)
```bash
conda create -n nbody python=3.10
conda activate nbody
pip install -r requirements.txt
conda install -c conda-forge cmake ninja
```
#### Option B — pyenv + venv (Linux / macOS / WSL)
```bash
pyenv install 3.10.14
pyenv local 3.10.14
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install cmake ninja
```

### Building the Engine
<div style="text-align: justify"> Always from the repository root, with the same env activated. </div>

| Goal | Command |
|-----------|--------|
| **Clean CPU build** | `python build.py --clean` |
| **Clean CPU + CUDA** | `python build.py --clean --cuda` |
| **Incremental CPU** | `python build.py` |
| **Incremental CUDA** | `python build.py --cuda` |

## USING THE SOLVER

To launch the GUI, in the repository root, with the same environment activated, do:
```bash
python main.py
```