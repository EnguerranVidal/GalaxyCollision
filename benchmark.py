from __future__ import annotations

import statistics
import sys
import time
import matplotlib.pyplot as plt

import numpy as np

sys.path.insert(0, ".")
import engine


SIZES = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000]
REPEATS = 10
WARMUP = 3
THETA = 1.0
G = 1.0
NEWTON_MAX = 50000
SEED = 0
SAVE_PATH = "benchmark_gpu_cpu.png"


def makeParticles(n: int, device: str):
    rng = np.random.default_rng(SEED)
    positionsArray = rng.uniform(-1.0, 1.0, size=(n, 3)).astype(np.float32)
    masses = np.ones(n, dtype=np.float32).tolist()
    group = engine.ParticleGroup(n, device)
    positions = [engine.Vector3(float(x), float(y), float(z)) for x, y, z in positionsArray]
    velocities = [engine.Vector3(0.0, 0.0, 0.0) for _ in range(n)]
    group.setPositions(positions)
    group.setVelocities(velocities)
    group.setMasses(masses)
    return group


def timeAverage(fn):
    try:
        for _ in range(WARMUP):
            fn()
        times = []
        for _ in range(REPEATS):
            t0 = time.perf_counter()
            fn()
            times.append(time.perf_counter() - t0)
        return statistics.mean(times)
    except Exception as exc:
        print(f"    [skip] {exc}")
        return None


def runBenchmark():
    newtonCalculator = engine.NewtonCalculator(G)
    barnesCalculator = engine.BarnesHutCalculator(THETA, G)
    results = {"NB_PARTICLES": [], "NEWTON_CPU": [], "NEWTON_GPU": [], "BARNES_CPU": [], "BARNES_GPU": []}
    for n in SIZES:
        print(f"measuring N = {n:,} ...", flush=True)
        results["NB_PARTICLES"].append(n)
        # NEWTONIAN CALCULATOR
        newtonCpuTime = newtonGpuTime = None
        if n <= NEWTON_MAX:
            cpu = makeParticles(n, "cpu")
            gpu = makeParticles(n, "gpu")
            newtonCpuTime = timeAverage(lambda: newtonCalculator.computeAccelerations(cpu))
            newtonGpuTime = timeAverage(lambda: newtonCalculator.computeAccelerations(gpu))
        else:
            print(f"  Newton skipped (N > {NEWTON_MAX})")
        results["NEWTON_CPU"].append(newtonCpuTime)
        results["NEWTON_GPU"].append(newtonGpuTime)
        # BARNES-HUT CALCULATOR
        barnesCpuTime = makeParticles(n, "cpu")
        barnesGpuTime = makeParticles(n, "gpu")
        results["BARNES_CPU"].append(timeAverage(lambda: barnesCalculator.computeAccelerations(barnesCpuTime)))
        results["BARNES_GPU"].append(timeAverage(lambda: barnesCalculator.computeAccelerations(barnesGpuTime)))
    return results


def plot(results) -> None:
    nbParticles = np.array(results["NB_PARTICLES"], dtype=float)
    series = [("NEWTON_CPU", "Newton CPU",  "C0", "-"), ("NEWTON_GPU", "Newton GPU",  "C0", "--"), ("BARNES_CPU", "Barnes CPU",  "C1", "-"), ("BARNES_GPU", "Barnes GPU",  "C1", "--")]
    fig, ax = plt.subplots(figsize=(9, 6))
    for key, label, color, ls in series:
        y = np.array([t if t is not None else np.nan for t in results[key]], dtype=float)
        mask = np.isfinite(y)
        if not mask.any():
            continue
        ax.plot(nbParticles[mask], y[mask], marker="o", color=color, linestyle=ls, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of particles N")
    ax.set_ylabel("Mean wall time (s)")
    ax.set_title("Force calculation: CPU vs GPU")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    if SAVE_PATH:
        fig.savefig(SAVE_PATH, dpi=150)
        print(f"saved: {SAVE_PATH}")
    plt.show()


def main():
    print(f"engine : {getattr(engine, '__file__', '?')}")
    print(f"sizes  : {SIZES}")
    print(f"repeats: {REPEATS}  |  theta={THETA}  |  newton_max={NEWTON_MAX}")
    results = runBenchmark()
    plot(results)


if __name__ == "__main__":
    main()