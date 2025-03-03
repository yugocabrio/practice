import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import traceback

try:
    # Read the CSV file
    print("Reading CSV file...")
    df = pd.read_csv('gkrfold_benchmark_results.csv')
    print("CSV file content:")
    print(df)

    # The CSV file already has the total values for GKR proofs and verification times
    # No need to multiply by n

    print("DataFrame after calculations:")
    print(df)

    # Create a figure with 3 subplots side by side
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # Plot 1: Comparison of gkrfold_time_ms and gkr_prove_time_ms
    ax1.plot(df['n'], df['gkrfold_time_ms'], 'r-', label='gkrfold_time_ms')
    ax1.plot(df['n'], df['gkr_prove_time_ms'], 'b-', label='gkr_prove_time_ms')
    ax1.set_xlabel('Number of Circuits (n)')
    ax1.set_ylabel('Time (ms)')
    ax1.set_title('Proving Time Comparison')
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.legend()

    # Plot 2: Comparison of folded_proof_size_kb and gkr_proof_size_kb
    ax2.plot(df['n'], df['folded_proof_size_kb'], 'r-', label='folded_proof_size_kb')
    ax2.plot(df['n'], df['gkr_proof_size_kb'], 'b-', label='gkr_proof_size_kb')
    ax2.set_xlabel('Number of Circuits (n)')
    ax2.set_ylabel('Size (KB)')
    ax2.set_title('Proof Size Comparison')
    ax2.set_xscale('log', base=2)
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.legend()

    # Plot 3: Comparison of gkrfold_verify_time_ms and gkr_verify_time_ms
    ax3.plot(df['n'], df['gkrfold_verify_time_ms'], 'r-', label='gkrfold_verify_time_ms')
    ax3.plot(df['n'], df['gkr_verify_time_ms'], 'b-', label='gkr_verify_time_ms')
    ax3.set_xlabel('Number of Circuits (n)')
    ax3.set_ylabel('Time (ms)')
    ax3.set_title('Verification Time Comparison')
    ax3.set_xscale('log', base=2)
    ax3.set_yscale('log')
    ax3.grid(True)
    ax3.legend()

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    plt.savefig('benchmark_comparison.png', dpi=300)

    # Show the figure
    plt.show()

    print("Graphs saved to benchmark_comparison.png")
except Exception as e:
    print(f"Error: {e}")
    traceback.print_exc()
