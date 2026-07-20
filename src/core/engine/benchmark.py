import time
import sys
import os
from random import randint

from src.core.engine.simulators import NBodySimulator
from src.core.engine.parameters import SimulatorParameters


def runBenchmark(device: str, nbParticles: int = 5000, nbSteps: int = 800):
    print(f"\n{'=' * 60}")
    print(f"RUNNING BENCHMARK: {device} | {nbParticles:,} particles | {nbSteps} steps")
    print(f"{'=' * 60}")
    parameters = SimulatorParameters()
    parameters.device = device
    parameters.distributionType = "BASIC"
    parameters.basicDistributionParameters.nbParticles = nbParticles
    parameters.seed = randint(0, 2 ** 31 - 1)
    parameters.timeStep = 0.01
    parameters.gravitationalConstant = 1.0
    simulator = NBodySimulator(parameters)
    simulator.setParameters(parameters)
    times = []
    startTotal = time.perf_counter()
    print(f"Starting simulation on {device}...")
    try:
        print(simulator.particles)
        for step in range(nbSteps):
            stepStartTime = time.perf_counter()
            simulator.step()
            stepTimeDuration = (time.perf_counter() - stepStartTime) * 1000
            times.append(stepTimeDuration)
            if step % 100 == 0 and step > 0:
                print(f"  Step {step:4d}/{nbSteps} - {stepTimeDuration:6.2f} ms")
    except Exception as e:
        print(f"Error during benchmark: {e}")
        return None
    totalTimeDuration = time.perf_counter() - startTotal
    averageTime = sum(times) / len(times)
    result = {"device": device, "particles": nbParticles, "steps": nbSteps, "avg_ms_per_step": round(averageTime, 3), "total_seconds": round(totalTimeDuration, 2), "estimated_fps": round(1000 / averageTime, 1), "min_ms": round(min(times), 3), "max_ms": round(max(times), 3)}
    print(f"\nRESULTS for {device}:")
    print(f"   Average time per step : {result['avg_ms_per_step']:6.3f} ms")
    print(f"   Estimated FPS         : {result['estimated_fps']:6.1f}")
    print(f"   Total time            : {result['total_seconds']:6.2f} seconds")
    print(f"   Min / Max step time   : {result['min_ms']:.3f} / {result['max_ms']:.3f} ms")
    return result


if __name__ == "__main__":
    print("N-Body Simulation Benchmark")
    print("===========================")
    nbParticlesList = [1000, 5000, 10000]
    steps = 600
    results = []
    for nbParticles in nbParticlesList:
        cpuResult = runBenchmark("CPU", nbParticles, steps)
        if cpuResult:
            results.append(cpuResult)
        time.sleep(1)
        gpuResult = runBenchmark("GPU", nbParticles, steps)
        if gpuResult:
            results.append(gpuResult)
        time.sleep(2)
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for r in results:
        print(
            f"{r['device']:4} | {r['particles']:6,} particles → {r['avg_ms_per_step']:6.3f} ms/step | ~{r['estimated_fps']:5.1f} FPS")
