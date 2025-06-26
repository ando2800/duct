import os
import pandas as pd
import matplotlib.pyplot as plt

# ==== Set font to Times New Roman for all plots ==== #
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Times New Roman'

# ==== Plotting functions ==== #
def plot_I(df, outputdir):
    plt.figure()
    plt.plot(df["time_nondim"], df["I"], color='black')
    plt.xlim(0.0, max(df["time_nondim"]))
    plt.ylim(-1.0, 1.0)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outputdir, "I_tsd.pdf"), bbox_inches='tight', format='pdf')

def plot_Ssum(df, outputdir):
    plt.figure()
    plt.plot(df["time_nondim"], df["Ssum"], color='black')
    plt.xlim(0.0, max(df["time_nondim"]))
    plt.ylim(0.0, 0.02)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outputdir, "Ssum_tsd.pdf"), bbox_inches='tight', format='pdf')

def plot_trajectory(df, outputdir):
    plt.figure(figsize=(6, 4))
    plt.plot(df["Ssum"], df["I"], color='black')
    plt.xticks([0.0, 0.005, 0.01, 0.015, 0.02], fontsize=10)
    plt.ylim(-1.0, 1.0)
    plt.xlim(0.0, 0.02)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outputdir, "trajectory_IvsSsum.pdf"), bbox_inches='tight', format='pdf')

def main(h5file, outputdir, plots=("I", "S", "trajectory")):
    os.makedirs(outputdir, exist_ok=True)

    # Load data from HDF5 file nessary for plotting
    df = pd.DataFrame()
    if "I" in plots or "trajectory" in plots:
        df["I"] = pd.read_hdf(h5file, "/new_statistics/Indicator")
    if "S" in plots or "trajectory" in plots:
        df["Ssum"] = pd.read_hdf(h5file, "/new_statistics/Ssum")
    if "I" in plots or "S" in plots or "trajectory" in plots:
        df["time_nondim"] = pd.read_hdf(h5file, "/time_nondim")

    # Check if the required columns are present
    if "I" in plots:
        plot_I(df, outputdir)
    if "S" in plots:
        plot_Ssum(df, outputdir)
    if "trajectory" in plots:
        plot_trajectory(df, outputdir)
    
    print("[MESSAGE] Selected plots have been saved to:", outputdir)

if __name__ == "__main__":
    h5file = "./data/processed/statistics_tsd.h5"
    outputdir = "./results/"
    plots = ("I", "S", "trajectory")  # (ex)：plots=("I",) だけにすると I のみ
    main(h5file, outputdir, plots)
