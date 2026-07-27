import argparse
import os
import sys
import time
from random import randint

projectRoot = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if projectRoot not in sys.path:
    sys.path.insert(0, projectRoot)

from src.core.engine.parameters import SimulatorParameters
from src.core.engine.simulators import NBodySimulator


def synchronizeDevice(particles):
    if particles.device != "gpu":
        return
    import cupy as cp
    cp.cuda.Stream.null.synchronize()


def runBenchmark(device: str, nbParticles: int = 5000, nbSteps: int = 200, integratorType: str = "EulerExplicit"):
    print(f"\n{'=' * 60}")
    print(f"RUNNING BENCHMARK: {device} | {nbParticles:,} particles | {nbSteps} steps | {integratorType}")
    print(f"{'=' * 60}")

    parameters = SimulatorParameters()
    parameters.device = device
    parameters.distributionType = "BASIC"
    parameters.basicDistributionParameters.nbParticles = nbParticles
    parameters.integratorType = integratorType
    parameters.calculatorType = "NEWTON"
    parameters.seed = randint(0, 2 ** 31 - 1)
    parameters.timeStep = 0.01
    parameters.gravitationalConstant = 1.0

    times = []
    print(f"Starting simulation on {device}...")
    try:
        simulator = NBodySimulator(parameters)
        if simulator.particles is None:
            raise RuntimeError("Simulator did not initialize particles")

        for _ in range(2):
            simulator.step()
        synchronizeDevice(simulator.particles)

        startTotal = time.perf_counter()
        for step in range(nbSteps):
            stepStartTime = time.perf_counter()
            simulator.step()
            synchronizeDevice(simulator.particles)
            stepTimeDuration = (time.perf_counter() - stepStartTime) * 1000
            times.append(stepTimeDuration)
            if step % 100 == 0 and step > 0:
                print(f"  Step {step:4d}/{nbSteps} - {stepTimeDuration:6.2f} ms")
    except Exception as e:
        print(f"Error during benchmark: {e}")
        return None

    totalTimeDuration = time.perf_counter() - startTotal
    averageTime = sum(times) / len(times)
    result = {
        "device": device,
        "particles": nbParticles,
        "steps": nbSteps,
        "avg_ms_per_step": round(averageTime, 3),
        "total_seconds": round(totalTimeDuration, 2),
        "estimated_fps": round(1000 / averageTime, 1),
        "min_ms": round(min(times), 3),
        "max_ms": round(max(times), 3),
    }
    print(f"\nRESULTS for {device}:")
    print(f"   Average time per step : {result['avg_ms_per_step']:6.3f} ms")
    print(f"   Estimated FPS         : {result['estimated_fps']:6.1f}")
    print(f"   Total time            : {result['total_seconds']:6.2f} seconds")
    print(f"   Min / Max step time   : {result['min_ms']:.3f} / {result['max_ms']:.3f} ms")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark the Galaxy Collision N-body engine.")
    parser.add_argument("--particles", nargs="+", type=int, default=[1000, 5000], help="Particle counts to benchmark.")
    parser.add_argument("--steps", type=int, default=200, help="Timed simulation steps per run.")
    parser.add_argument("--devices", nargs="+", choices=["CPU", "GPU", "cpu", "gpu"], default=["CPU", "GPU"], help="Devices to benchmark.")
    parser.add_argument("--integrator", choices=["EulerExplicit", "RK4"], default="EulerExplicit", help="Integrator used during the benchmark.")
    args = parser.parse_args()

    print("N-Body Simulation Benchmark")
    print("===========================")
    results = []
    for nbParticles in args.particles:
        for device in args.devices:
            result = runBenchmark(device.upper(), nbParticles, args.steps, args.integrator)
            if result:
                results.append(result)
            time.sleep(1)

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for result in results:
        print(
            f"{result['device']:4} | {result['particles']:6,} particles -> "
            f"{result['avg_ms_per_step']:6.3f} ms/step | ~{result['estimated_fps']:5.1f} FPS"
        )
