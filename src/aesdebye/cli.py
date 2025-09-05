import argparse
import sys

import ase
import ase.io
import matplotlib.pyplot as plt

import aesdebye

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        prog="aesdebye", description="Debye implementation CLI"
    )

    # --- INPUT CONFIGURATION ---
    input_group = parser.add_argument_group("Input configuration")
    input_group.add_argument(
        "-f",
        "--inputFilename",
        type=str,
        default="",
        help="Input filename, any format supported by the ase.io.read function",
    )

    # --- OUTPUT OPTIONS ---
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "-o",
        "--outputDir",
        type=str,
        default="",
        help="Output directory to save the data",
    )
    output_group.add_argument(
        "-p", "--plot", action="store_true", default=False, help="Plot the results"
    )
    output_group.add_argument(
        "-nv",
        "--nonVerbose",
        action="store_true",
        default=False,
        help="Disable verbose output",
    )

    # --- CPU OPTIONS ---
    cpu_group = parser.add_argument_group("CPU options")
    cpu_group.add_argument(
        "-nt", "--nThreads", type=int, default=-1, help="Number of threads to use"
    )
    cpu_group.add_argument(
        "-nc", "--nCells", type=int, default=15, help="Number of cells in the lattice"
    )
    cpu_group.add_argument(
        "-nlh",
        "--dontUseLocalHistogram",
        action="store_true",
        default=False,
        help="Dont use local histogram for noisy calcs",
    )

    # --- GPU OPTIONS ---
    gpu_group = parser.add_argument_group("GPU options")
    gpu_group.add_argument(
        "-gpu",
        "--useGPU",
        action="store_true",
        default=False,
        help="Use GPU for parallelization",
    )
    gpu_group.add_argument(
        "-pc",
        "--pseudoCoal",
        action="store_true",
        default=False,
        help="Use pseudo coal",
    )
    gpu_group.add_argument(
        "-fg",
        "--fillGPU",
        action="store_true",
        default=False,
        help="Fill the GPU with threads",
    )

    # --- MULTI-CPU/GPU OPTIONS ---
    common_group = parser.add_argument_group("Common options: CPU/GPU")
    common_group.add_argument(
        "-mpi",
        "--useMPI",
        action="store_true",
        default=False,
        help="Use MPI for parallelization",
    )
    common_group.add_argument(
        "-br",
        "--binsResolution",
        type=int,
        default=1.0,
        help="Number of bins in the histogram",
    )

    # --- RANGE & PHYSICS SETTINGS ---
    physics_group = parser.add_argument_group("Range and physics settings")
    physics_group.add_argument(
        "-st", "--start", type=float, default=0.0, help="Start of the q/theta range"
    )
    physics_group.add_argument(
        "-e", "--end", type=float, default=20.0, help="End of the q/theta range"
    )
    physics_group.add_argument(
        "-stps",
        "--steps",
        type=int,
        default=2000,
        help="Number of steps in the q/theta range",
    )
    physics_group.add_argument(
        "-wl",
        "--wavelength",
        type=float,
        default=0.4,
        help="Wavelength. This will be used to convert between q and two theta. q = 4pi/wavelength * sin(theta)",
    )
    physics_group.add_argument(
        "-tt",
        "--twoThetaSpace",
        action="store_true",
        default=False,
        help="Use two theta instead of q",
    )

    # --- BENCHMARK CONFIGURATION ---
    benchmark_group = parser.add_argument_group("Benchmark configuration")
    benchmark_group.add_argument(
        "-nr",
        "--nRepeats",
        type=int,
        default=-1,
        help="Number of repeats in the lattice",
    )
    benchmark_group.add_argument(
        "-s",
        "--stdDev",
        type=float,
        default=0.0,
        help="Std dev of noise in the lattice",
    )

    return parser.parse_args(), parser


def main():
    args, parser = parse_args()

    # Ensure nRepeats > 0 or inputFilename is provided
    if args.nRepeats < 0 and not args.inputFilename:
        print("Either nRepeats or inputFilename must be provided", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # Initialize DebyeCalculator
    calc = aesdebye.DebyeCalculator(
        nThreads=args.nThreads,
        nCells=args.nCells,
        binsResolution=args.binsResolution,
        useMPI=args.useMPI,
        useGPU=args.useGPU,
        useLocalHist=not args.dontUseLocalHistogram,
        smallBin=False,
        pseudoCoal=args.pseudoCoal,
        fillGPU=args.fillGPU,
        verbose=not args.nonVerbose,
    )

    # Load positions
    if args.nRepeats > 0:
        positions = aesdebye.generateData(3.89070, args.nRepeats, args.stdDev, 10)
    else:
        atoms = ase.io.read(args.inputFilename)
        positions = aesdebye.Positions(
            chemicalSymbols=atoms.get_chemical_symbols(),
            coordinates=atoms.get_positions(),
        )

    results = calc.calculateProfile(
        positions,
        args.start,
        args.end,
        args.steps,
        args.twoThetaSpace,
        args.wavelength,
        "",
    )

    if not args.nonVerbose and calc.parallelHelper.world_rank == 0:
        last_pdf, last_profile = list(results.values())[-1]
        print(last_profile)

    if args.plot and calc.parallelHelper.world_rank == 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        for name, (_, profile) in results.items():
            ax.semilogy(profile.q, profile.intensity, label=name)

        ax.set_xlabel(
            r"$q$ [$\AA^{-1}$]" if not args.twoThetaSpace else r"$2\theta$ [deg]"
        )
        ax.set_ylabel("Intensity [a.u.]")
        ax.legend()
        plt.show()

    if args.outputDir and calc.parallelHelper.world_rank == 0:
        prefix = (
            args.inputFilename.rsplit(".", 1)[0]
            if args.inputFilename
            else f"{args.nRepeats}_repeats_{args.stdDev}_stdDev_"
        )
        for name, (pdf, profile) in results.items():
            output_prefix = f"{args.outputDir}/{prefix}_{name}_"
            pdf.toCSV(output_prefix + "_pdf.csv")
            profile.toCSV(output_prefix + "_profile.csv")


if __name__ == "__main__":
    main()
