import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_pfm(file):
    with open(file, 'rb') as f:
        header = f.readline().decode('utf-8').rstrip()
        if header not in ['PF', 'Pf']:
            raise Exception('Not a valid PFM file.')
        color = (header == 'PF')

        dims = f.readline().decode('utf-8')
        while dims.startswith('#'):
            dims = f.readline().decode('utf-8')
        width, height = map(int, dims.strip().split())

        scale = float(f.readline().decode('utf-8').strip())
        endian = '<' if scale < 0 else '>'
        
        data = np.fromfile(f, endian + 'f')
        shape = (height, width, 3) if color else (height, width)
        data = np.reshape(data, shape)
        return np.flipud(data)

def compute_rmse(img, ref):
    if img.shape != ref.shape:
        raise ValueError("Image and reference must have the same shape.")
    return np.sqrt(np.mean((img - ref) ** 2))

def extract_spp(filename):
    match = re.search(r'seq_(\d+)', filename)
    return int(match.group(1)) if match else None

def extract_m_multiplier(folder_name):
    """
    Extracts the number after 'M' in folder names like 'C8-M8' or 'Hero-M16'.
    """
    match = re.search(r'M(\d+)', folder_name)
    if match:
        return int(match.group(1))
    return 1  # Default to 1 if M is not found

def main(root_path, reference_name="reference.pfm"):
    reference_path = os.path.join(root_path, reference_name)
    if not os.path.exists(reference_path):
        print(f"Error: Reference image not found: {reference_path}")
        sys.exit(1)

    print(f"Loading reference: {reference_name}")
    reference = read_pfm(reference_path)

    plt.figure(figsize=(10, 6))
    
    methods = [d for d in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, d))]
    
    for method in sorted(methods):
        method_path = os.path.join(root_path, method)
        spp_accumulator = {}
        
        # Get multiplier from folder name (e.g., M8 -> 8)
        m_multiplier = extract_m_multiplier(method)

        runs = [d for d in os.listdir(method_path) if d.startswith("run-") and os.path.isdir(os.path.join(method_path, d))]
        if not runs:
            continue

        print(f"Processing method: {method} (M={m_multiplier}, {len(runs)} runs)")

        for run in runs:
            run_path = os.path.join(method_path, run)
            for file in os.listdir(run_path):
                if file.endswith(".pfm"):
                    spp = extract_spp(file)
                    if spp is not None:
                        try:
                            img = read_pfm(os.path.join(run_path, file))
                            rmse = compute_rmse(img, reference)
                            
                            if spp not in spp_accumulator:
                                spp_accumulator[spp] = []
                            spp_accumulator[spp].append(rmse)
                        except Exception as e:
                            print(f"  Error reading {file}: {e}")

        if spp_accumulator:
            # Sort original SPP values
            sorted_spps = np.array(sorted(spp_accumulator.keys()))
            avg_rmses = [np.mean(spp_accumulator[s]) for s in sorted_spps]
            
            # Multiply SPP by the M factor for the final plot x-axis
            total_samples = sorted_spps * m_multiplier
            
            plt.loglog(total_samples, avg_rmses, marker='o', label=method)

    plt.xlabel("Total Samples (SPP * M)")
    plt.ylabel("Average RMSE")
    plt.title("RMSE Convergence vs Total Samples (log-log scale)")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    output_file = os.path.join(root_path, "averaged_rmse_plot.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Plot saved to: {output_file}")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 plot_avg_rmse.py <path_to_hero_folder> [reference_filename]")
        sys.exit(1)

    path = sys.argv[1]
    ref_file = sys.argv[2] if len(sys.argv) > 2 else "reference.pfm"
    main(path, ref_file)