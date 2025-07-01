#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_probability(path):
    data = np.loadtxt(path)
    x = data[:, 0].astype(int)
    y = data[:, 1].astype(int)
    prob = data[:, 4]
    nx = x.max() + 1
    ny = y.max() + 1
    grid = np.zeros((ny, nx))
    for xi, yi, p in zip(x, y, prob):
        grid[yi, xi] = p
    return grid


def main():
    parser = argparse.ArgumentParser(
        description="Plot probability density from a slice_fields data file"
    )
    parser.add_argument("input", help="Path to slices_fields_#.dat file")
    parser.add_argument(
        "--output",
        "-o",
        help="Optional path to save the plot as PNG. If omitted, show window.",
    )
    args = parser.parse_args()

    grid = load_probability(args.input)
    plt.imshow(grid, cmap="viridis", origin="lower")
    plt.colorbar(label="Probability")
    plt.xlabel("X index")
    plt.ylabel("Y index")
    plt.title("Probability density")

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()


if __name__ == "__main__":
    main()
