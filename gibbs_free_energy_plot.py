import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker

# Set publication-quality plot parameters
plt.rcParams.update({
    "text.usetex": False,               # Use MathText (not external LaTeX)
    "svg.fonttype": "none",             # Keep text editable in SVG
    "font.family": "Times New Roman",   # Use Times New Roman
    "font.size": 16,
    "mathtext.fontset": "stix",            # Optional: math style (cm, stix, etc.)
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "lines.linewidth": 2.0,
    "axes.linewidth": 1.5,
})





def plot_energy_profile(ax, delta_G, color="b", linestyle="--", label=None, show_dG=True):
    """
    Plot one Gibbs free energy profile on a given axis.
    
    Parameters:
    -----------
    ax : matplotlib axis
        Axis to plot on.
    delta_G : list of floats
        Gibbs free energy changes for each step.
    intermediates : list of str
        Names of intermediates (length must be len(delta_G)+1).
    color : str
        Color for horizontal energy levels.
    linestyle : str
        Style for connecting lines.
    label : str
        Label for legend.
    """
    # Compute cumulative Gibbs energies
    G_levels = [0]  # reference starts at 0
    for dG in delta_G:
        G_levels.append(G_levels[-1] + dG)

    # Draw energy levels and connecting lines
    for i, G in enumerate(G_levels):
        ax.hlines(G, i-0.2, i+0.2, colors=color, linewidth=3)
        if i < len(G_levels)-1:
            ax.plot([i+0.2, i+1-0.2], [G, G_levels[i+1]], color=color, linestyle=linestyle)
            if show_dG:
                mid_x = i + 0.5
                mid_y = (G + G_levels[i+1]) / 2
                ax.text(mid_x, mid_y + 0.1, f"{delta_G[i]:+.2f}", 
                        ha="center", va="bottom", fontsize=16, color=color)
    
    if label:
        ax.plot([], [], color=color, linewidth=3, label=label)

    return G_levels


# ---------------- Example usage ---------------- #

intermediates = [r'$*Na_2S_8$', r'$*Na_2S_6$', r'$*Na_2S_4$', r'$*Na_2S_2$', r'$*Na_2S$']

num_intermediates = len(intermediates)
X_levels = [r'$*S_8$'] + intermediates

delta_G1 = [-0.5, 0.8, -0.3, 0.6, -0.4]
delta_G2 = [-0.2, 0.5, -0.1, 0.3, -0.6]  # another pathway
with open("energy_results.dat") as energy_results:
    lines = energy_results.readlines()
    for line in lines:
        if "gibbs free energy list:" in line:
            print(line)
            delta_Gs_str = line.split("gibbs free energy list:")[1]
            delta_Gs = list(map(float, delta_Gs_str.split()))




fig, ax = plt.subplots(figsize=(9,3))

# Plot multiple profiles
surfaces = ["g-CN", "m-CN", "g-C", "N3g-C", "N4g-C"]

# Plot multiple profiles
colors = ["r", "g", "b", "orange", "m"]

for i in range(len(surfaces)):
    print(i)
    delta_G = delta_Gs[i*num_intermediates : num_intermediates*(i+1)]
    print(delta_G)
    plot_energy_profile(ax, delta_G, color=colors[i], linestyle="--", label="{}".format(surfaces[i]))


# Common formatting
ax.set_xticks(range(len(X_levels)))
ax.set_xticklabels(X_levels, ha="center")
ax.tick_params(axis="x", which="both", length=0)

ax.set_ylabel("Gibbs Free Energy (ΔG)")
ax.set_title("Na-S battery discharging curve")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)

ax.legend()

plt.tight_layout()
plt.savefig("reaction_energy_profile.svg", format="svg")
plt.show()
