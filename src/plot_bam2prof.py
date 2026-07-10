import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

parser = argparse.ArgumentParser(description='Plot 5\' and 3\' substitution rate profiles produced by bam2prof.')
parser.add_argument('fivep_end_file_path', help='Path to the 5\' end profile file')
parser.add_argument('threep_end_file_path', help='Path to the 3\' end profile file')
parser.add_argument('output_prefix', nargs='?', default='bam2prof', help='Prefix for the output plot files (default: bam2prof)')
parser.add_argument('--title', default=None, help='Custom title to use for the plots')
parser.add_argument('--xlim', type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=None, help='Fixed x-axis range instead of auto-detecting it')
parser.add_argument('--ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'), default=None, help='Fixed y-axis range instead of auto-detecting it')
args = parser.parse_args()

fivep_end_file_path  = args.fivep_end_file_path
threep_end_file_path = args.threep_end_file_path
output_prefix = args.output_prefix
custom_title = args.title
custom_xlim = args.xlim
custom_ylim = args.ylim

# Read the data
fivep_data = pd.read_csv(fivep_end_file_path, delimiter='\t')
threep_data = pd.read_csv(threep_end_file_path, delimiter='\t')

# Reverse the order of the threep_data to have the first line as the leftmost position
threep_data = threep_data.iloc[::-1].reset_index(drop=True)

# Define line widths and colors
thick_line_width = 2.5  # Thicker line for specific substitutions
thin_line_width = 1.0   # Thinner line for the rest


# Function to plot data with adjustments for y-axis scale, label position, improved legend placement, and larger legend
def plot_data(data, title, reverse_xaxis=False, save_path='', ylim=None, xlim=None, yticks_right=False):
    fig, ax = plt.subplots(figsize=(10, 8))  # Create a new figure for each plot
    lines, labels = [], []  # Track labels for the legend
    for column in data.columns:
        if column == 'C>T':
            line, = ax.plot(data.index, data[column], label=column, color='red', linewidth=thick_line_width)
        elif column == 'G>A':
            line, = ax.plot(data.index, data[column], label=column, color='blue', linewidth=thick_line_width)
        else:
            line, = ax.plot(data.index, data[column], label=column, linewidth=thin_line_width)
        lines.append(line)
    labels = [l.get_label() for l in lines]  # Update labels from the last plot
    ax.set_xlabel('Position on the fragment (bp)')
    ax.set_title(title)
    ax.grid(True)
    if reverse_xaxis:
        ax.set_xticks(np.arange(len(data)))
        ax.set_xticklabels(np.arange(len(data))[::-1], rotation=90)
    else:
        ax.set_xticks(np.arange(len(data)))
        ax.set_xticklabels(np.arange(len(data)), rotation=90)
    if ylim:
        ax.set_ylim(ylim)  # Apply the determined y-axis limits
    if xlim:
        ax.set_xlim(xlim)  # Apply the fixed x-axis range
    if yticks_right:
        ax.yaxis.tick_right()  # Move y-axis labels to the right
        ax.yaxis.set_label_position("right")
    # Adjust legend placement and size
    ax.legend(lines, labels, loc='center left' if yticks_right else 'upper left',  #yticks_right=true if 3'
              bbox_to_anchor=(1.1, 0.5) if yticks_right else (1.05, 0.7),  #yticks_right=true if 3'
              borderaxespad=0., fontsize='large')  # Increase the legend font size
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(save_path)
    plt.close(fig)

# Determine the common y-axis scale, unless the user supplied a fixed range
all_data = pd.concat([fivep_data, threep_data], axis=0)
common_ylim = custom_ylim if custom_ylim else (all_data.min().min(), min(all_data.max().max()+0.05,1.0))

fivep_title = f'{custom_title} - 5\' Substitution Rates' if custom_title else '5\' Substitution Rates'
threep_title = f'{custom_title} - 3\' Substitution Rates' if custom_title else '3\' Substitution Rates'

# Plot and save left end data with common y-axis scale
plot_data(fivep_data, fivep_title, save_path=f'{output_prefix}_5p_end_plot.pdf', ylim=common_ylim, xlim=custom_xlim)

# Plot and save right end data with common y-axis scale and y-axis labels on the right
plot_data(threep_data, threep_title, reverse_xaxis=True, save_path=f'{output_prefix}_3p_end_plot.pdf', ylim=common_ylim, xlim=custom_xlim, yticks_right=True)


# Create a single figure with two subplots side by side
fig, axs = plt.subplots(1, 2, figsize=(20, 8))


# Function to plot data on given axis
def plot_data_combined(ax, data, title, reverse_xaxis=False, ylim=None, xlim=None, yticks_right=False):
    lines, labels = [], []  # Track labels for the legend
    for column in data.columns:
        if column == 'C>T':
            line, = ax.plot(data.index, data[column], label=column, color='red', linewidth=thick_line_width)
        elif column == 'G>A':
            line, = ax.plot(data.index, data[column], label=column, color='blue', linewidth=thick_line_width)
        else:
            line, = ax.plot(data.index, data[column], label=column, linewidth=thin_line_width)
        lines.append(line)
    ax.set_xlabel('Position on the fragment (bp)')
    ax.set_title(title)
    ax.grid(True)
    if reverse_xaxis:
        ax.set_xticks(np.arange(len(data)))
        ax.set_xticklabels(np.arange(len(data))[::-1], rotation=90)
    else:
        ax.set_xticks(np.arange(len(data)))
        ax.set_xticklabels(np.arange(len(data)), rotation=90)
    if ylim:
        ax.set_ylim(ylim)
    if xlim:
        ax.set_xlim(xlim)
    if yticks_right:
        ax.yaxis.tick_right()
    return lines, labels



# Define line widths and colors as before
thick_line_width = 2.5  # Thicker line for specific substitutions
thin_line_width  = 1.0   # Thinner line for the rest

# Determine the common y-axis scale, unless the user supplied a fixed range
# Assuming all_data is a combination of fivep_data and right_data
common_ylim = custom_ylim if custom_ylim else (all_data.min().min(), min(all_data.max().max()+0.05,1.0))

# Adjust the layout to include a third axis for the legend at the top or bottom
fig, axs = plt.subplots(1, 3, figsize=(20, 8), gridspec_kw={'width_ratios': [1, 1, 0.05]})

# Hide the third axis (used for legend) by making it invisible
axs[2].axis('off')  # This axis will not display any plot, just the legend

# Plot data on the fivep and right axis
lines1, labels1 = plot_data_combined(axs[0], fivep_data, fivep_title,                     ylim=common_ylim, xlim=custom_xlim)
lines2, labels2 = plot_data_combined(axs[1], threep_data, threep_title, reverse_xaxis=True, ylim=common_ylim, xlim=custom_xlim, yticks_right=True)

# Since the lines and labels are the same for both plots, use one set for the legend
# Place the legend in the invisible third axis
handles, labels = axs[0].get_legend_handles_labels()  # Assuming same labels for both plots
axs[2].legend(handles, labels, loc='center left', ncol=1, borderaxespad=0.)


plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to make space for the legend
plt.subplots_adjust(right=0.95)  # Adjust this value as needed to reduce white space

plt.savefig(f'{output_prefix}_combined_plot_with_legend.pdf')
plt.close(fig)
