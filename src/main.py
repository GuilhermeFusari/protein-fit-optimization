#!/usr/bin/env python3
"""
Main entry point for the PDB alignment and packing pipeline.
Handles argument parsing and routes to the appropriate mode.
"""
import argparse
import sys
import traceback
from multiprocessing import cpu_count

try:
    from icp_aligner import find_best_pdb
    from pdb_packing import run_packing_pipeline
except ImportError as e:
    print(f"Import error: {e}")
    print("Ensure you run from the project root: python3 src/main.py ...")
    sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(description='PDB Alignment and Packing Tool')
    parser.add_argument('mode', choices=['single', 'packing'], help="Execution mode.")
    parser.add_argument('-i', '--input', required=True, help='Input folder or single PDB file.')
    parser.add_argument('-e', '--envelope', required=True, help='Target CIF envelope file.')
    parser.add_argument('-o', '--output', required=True, help='Output directory.')
    
    # Algorithmic parameters
    parser.add_argument('--max-iter', type=int, default=50, help='Max ICP iterations.')
    parser.add_argument('--max-structures', type=int, default=20, help='Max structures for packing mode.')
    parser.add_argument('--penalty', type=float, default=50.0, help='Boundary violation penalty weight (packing only).')
    
    # Performance & System args
    parser.add_argument('--workers', type=int, default=cpu_count(), help='Number of CPU threads.')
    parser.add_argument('--sample-env', type=int, default=5000, help='Envelope point downsampling limit.')
    parser.add_argument('--align-what', default='protein', choices=['protein', 'envelope'], help='Alignment target.')

    return parser.parse_args()

def main():
    args = parse_arguments()
    print(f"Initializing mode: {args.mode.upper()}")

    if args.mode == 'single':
        try:
            find_best_pdb(args)
        except Exception as e:
            print(f"Failed in single mode: {e}")
            traceback.print_exc()
            sys.exit(1)
            
    elif args.mode == 'packing':
        try:
            run_packing_pipeline(args.input, args.envelope, args.output, args)
        except Exception as e:
            print(f"Failed in packing mode: {e}")
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()