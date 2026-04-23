#!/usr/bin/env python3

import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="+", required=True)
    parser.add_argument("--iqtrees", nargs="+", required=True)
    parser.add_argument("--genes", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def extract_stats(iqtree_file, best_model):
    logl = aic = bic = "NA"
    with open(iqtree_file) as f:
        for line in f:
            # Try labeled lines first (some IQ-TREE versions use these)
            if "Log-likelihood of the tree:" in line:
                logl = line.strip().split()[-1]
            elif "Akaike information criterion (AIC) score:" in line:
                aic = line.strip().split()[-1]
            elif "Bayesian information criterion (BIC) score:" in line:
                bic = line.strip().split()[-1]

            # Parse from model selection table row matching best-fit model
            # Table format: No. Model -LnL df AIC AICc BIC
            stripped = line.strip()
            if best_model and best_model in stripped:
                parts = stripped.split()
                # Expected: index, model_name, LnL, df, AIC, AICc, BIC
                if len(parts) >= 7:
                    try:
                        logl = str(round(float(parts[2]), 3))
                        aic  = str(round(float(parts[4]), 3))
                        bic  = str(round(float(parts[6]), 3))
                    except ValueError:
                        pass

    return logl, aic, bic


def main():
    args = parse_args()

    if not (len(args.models) == len(args.iqtrees) == len(args.genes)):
        raise ValueError("Mismatch: models, iqtrees, and genes must have same length")

    with open(args.output, "w") as out:
        out.write("Best-fit Substitution Models Summary\n")
        out.write("=====================================\n\n")
        out.write("Gene\tBest_Model\tLogL\tAIC\tBIC\n")

        for gene, model_file, iqtree_file in zip(
            args.genes, args.models, args.iqtrees
        ):
            with open(model_file) as mf:
                model = mf.read().strip()

            logl, aic, bic = extract_stats(iqtree_file, model)

            out.write(f"{gene}\t{model}\t{logl}\t{aic}\t{bic}\n")


if __name__ == "__main__":
    main()
