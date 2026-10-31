#!/usr/bin/env python3
# USAGE:
# ./script min_comp max_con infile1.stats [file2.stats file3.stats]

import sys
import os
import matplotlib.pyplot as plt

plt.switch_backend('agg')

min_completion = int(sys.argv[1])
max_contamination = int(sys.argv[2])

file_list = sys.argv[3:]


####################################################################################################################################
############################################         MAKE THE COMPLETION PLOT           ############################################
####################################################################################################################################
print("Loading completion info....")
data = {}
max_x = 0
# loop over all bin .stats files
for file_name in file_list:
    print(f"Loading {file_name}...")
    bin_set = ".".join(file_name.split("/")[-1].split(".")[:-1])
    data[bin_set] = []
    try:
        with open(file_name) as f:
            header = next(f)
            seperator_line = next(f)
            for line in f:
                bin_data = line.strip().split()

                try:
                    bin_comp = float(bin_data[1])
                    bin_cont = float(bin_data[2])
                except IndexError:
                    print(f"Error: {line}, file: {file_name}")
                    raise IndexError
                
                # skip bins that are too contaminated or very incomplete
                try:
                    if bin_comp < min_completion:
                        continue
                    if bin_cont > max_contamination:
                        continue
                except IndexError:
                    print(f"Error: {line}, file: {file_name}")
                    raise IndexError

                # save the completion value of each bin into a list
                data[bin_set].append(bin_comp)
        
        if len(data[bin_set]) > max_x:
            max_x = len(data[bin_set])
    except StopIteration:
        print(f"Warning: {file_name} is empty. Skipping...")
        continue

# Pass if no bins are found
if all([len(data[k]) == 0 for k in data]):
    print("No bins found. Exiting...")
    sys.exit(0)

# sort the completion data sets
for bin_set in data:
    data[bin_set].sort(reverse=True)

print("Plotting completion data...")
# MAKING THE PLOT PRETTY!!!!
# set some color schemes
tableau20 = [(214, 39, 40), (31, 119, 180), (255, 127, 14),
             (44, 160, 44), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

tableau20 = [(r / 255., g / 255., b / 255.) for r, g, b in tableau20]
plot_colors = {}
for i, label in enumerate(file_list):
    bin_set = ".".join(label.split("/")[-1].split(".")[:-1])
    plot_colors[bin_set] = tableau20[i % len(tableau20)]

# set figure size
plt.figure(figsize=(16, 8))
plt.style.use('ggplot')

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(121)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines['bottom'].set_color('black')
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
plt.ylim(min_completion, 105)
max_x = 0
for k in data:
    if len(data[k]) > max_x: max_x = len(data[k])
plt.xlim(0, max_x)

# Make sure your axis ticks are large enough to be easily read.
plt.yticks(range(min_completion, 105, 10), [str(x) + "%" for x in range(min_completion, 105, 10)], fontsize=14)
plt.xticks(fontsize=14)

# Provide tick lines across the plot to help your viewers trace along
for y in range(min_completion, 105, 10):
    plt.axhline(y=y, linestyle="--", lw=0.5, color="black", alpha=0.3)
for x in range(0, 1000, 20):
    plt.axvline(x=x, linestyle="--", lw=0.5, color="black", alpha=0.3)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
plt.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=True, left=False, right=False, labelleft=True)

# PLOTTING THE DATA
# prepare labels
labels = []
for k in data: labels.append(k)

# plot the data and labels
N = len(labels)
y_increment = (100 - min_completion) / N / 2
y_pos = 100 - y_increment

for rank, bin_set in enumerate(labels):
    # chose a color!
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot
    for x_pos, y in enumerate(data[bin_set]):
        if y < y_pos:
            break
        plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
        y_pos -= y_increment

# add plot and axis titles and adjust edges
plt.title("Bin completion ranking", fontsize=26)
plt.xlabel("Descending completion rank", fontsize=16)
plt.ylabel("Estimated bin completion", fontsize=16)

####################################################################################################################################
############################################         MAKE THE CONTAMINATION PLOT        ############################################
####################################################################################################################################
print("Loading contamination info...")

data = {}
# loop over all bin .stats files
for file_name in file_list:
    print(f"Loading {file_name}...")
    bin_set = ".".join(file_name.split("/")[-1].split(".")[:-1])
    data[bin_set] = []
    with open(file_name) as f:
        header = next(f)
        seperator_line = next(f)
        for line in f:
            bin_data = line.strip().split()
            try:
                bin_comp = float(bin_data[1])
                bin_cont = float(bin_data[2])
            except IndexError:
                print(f"Error: {line}, file: {file_name}")
                raise IndexError
            # skip bins that are too incomplete or way too contaminated
            try:
                if bin_comp < min_completion:
                    continue
                if bin_cont > max_contamination:
                    continue
            except IndexError:
                print(f"Error: {line}, file: {file_name}")
                raise IndexError

            # save the contamination value of each bin into a list
            data[bin_set].append(bin_cont)

# sort the contamination data sets
for bin_set in data:
    data[bin_set].sort(reverse=False)

print("Plotting the contamination data...")
# MAKING THE PLOT PRETTY!!!!
# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(122)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines['bottom'].set_color('black')
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
plt.ylim(0, max_contamination + 1)
max_x = 0
for k in data:
    if len(data[k]) > max_x:
        max_x = len(data[k])
if max_x == 0:
    print(f"data: {data}")
    max_x = 1
plt.xlim(0, max_x)

# Make sure your axis ticks are large enough to be easily read.
plt.yticks(range(0, max_contamination + 1, 1), [str(x) + "%" for x in range(0, max_contamination + 1, 1)], fontsize=14)
plt.xticks(fontsize=14)

# Provide tick lines across the plot to help your viewers trace along
for y in range(0, max_contamination + 1, 1):
    plt.axhline(y=y, linestyle="--", lw=0.5, color="black", alpha=0.3)
for x in range(0, 1000, 20):
    plt.axvline(x=x, linestyle="--", lw=0.5, color="black", alpha=0.3)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
plt.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=True, left=False, right=False, labelleft=True)

# PLOTTING THE DATA
# prepare labels
labels = []
for k in data: labels.append(k)

# plot the data and labels
N = len(labels)
y_increment = max_contamination / N / 2
y_pos = y_increment

for rank, bin_set in enumerate(labels):
    # chose a color!
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot
    for x_pos, y in enumerate(data[bin_set]):
        if y > y_pos:
            break
        plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
        y_pos += y_increment

# add plot and axis titles and adjust the edges
plt.title("Bin contamination ranking", fontsize=26)
plt.xlabel("Ascending contamination rank", fontsize=16)
plt.ylabel("Estimated bin contamination (log scale)", fontsize=16)
plt.gcf().subplots_adjust(right=0.9)

# save figure
print("Saving figure binning_results.png ...")
plt.tight_layout(w_pad=10)
plt.subplots_adjust(top=0.92, right=0.90, left=0.08)
plt.savefig("binning_results.png", format='png', dpi=600)
