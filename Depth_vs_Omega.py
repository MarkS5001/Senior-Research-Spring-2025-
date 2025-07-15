import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

def plot_omega_vs_depth(path_name: str):

    data = pd.read_excel(path_name)

    omega = data['omega']
    depth_curve_fit = data['effective_depth']
    depth_intensity = data['effective_depth_intensity']
    depth_summation = data['effective_depth_sum']

    fig, ax = plt.subplots()
    ax.plot(omega, depth_curve_fit, label='Curve Fit', marker='o')
    ax.plot(omega, depth_intensity, label='Intensity', marker='o')
    ax.plot(omega, depth_summation, label='Summation', marker='o')
    ax.legend()
    ax.set_title('Omega vs Depth')
    ax.set_xlabel('Omega (°)')
    ax.set_ylabel(' Depth (μm)')
    ax.set_xticks(omega[5:], rotation=90)
    ax.set_xticklabels(omega[5:], rotation=90)
    ax.set_yticks(depth_curve_fit[5::2])
    ax.invert_yaxis() # Flip the y-axis so that depth increases downwards

    # Add inset axes in the lower left corner (loc=3)
    ax_inset = inset_axes(ax, width="40%", height="30%", loc=3, bbox_to_anchor=(0.1, 0.13, 1, 1), bbox_transform=ax.transAxes)
    ax_inset.plot(omega, depth_curve_fit, marker='o')
    ax_inset.plot(omega, depth_intensity, marker='o')
    ax_inset.plot(omega, depth_summation, marker='o')
    # ax_inset.set_title("Small Omega")

    ax_inset.set_xlim(-0.025, .6)
    ax_inset.set_ylim(-0.2,4)
    ax_inset.set_xticks(omega[:5], rotation=90)
    ax_inset.set_xticklabels(omega[:5], rotation=90)
    ax_inset.set_yticks(depth_curve_fit[:5:2])
    ax_inset.invert_yaxis() # Flip the y-axis so that depth increases downwards

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    material = input("p, p01, or dt: ")

    path_name = f'E://Full//{material} csv' # input("Enter the path to the XRD data folder: ").strip()  # Get the path to the XRD data files from the user

    try:
        plot_omega_vs_depth(path_name + '//Processed_XRD_Data.xlsx')
    
    except FileNotFoundError as e:
        print(e)