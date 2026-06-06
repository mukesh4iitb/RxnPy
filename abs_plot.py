import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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




# Replace this with the path to your file


if os.path.exists("abs_plot.txt"):
    df = pd.read_csv("abs_plot.txt", sep="\s+", comment="#")
else:
    print("\nCreate abs_plot.txt file to plot absorption energy with legends.\n")
    with open("energy_results.dat") as energy_results:
        lines = energy_results.readlines()
        for line in lines:
            if "Absorption energy list:" in line:
                print(line)
                abs_en_str = line.split("Absorption energy list:")[1].split()
                #print(abs_en_str)
                intermediate = 5 #input("Enter number of intermediate:\n")
                intermediate_sys_list =[abs_en_str[i:i+intermediate] for i in range(0, len(abs_en_str), intermediate)]  # intermediate as row and sys as column.
                df=pd.DataFrame(intermediate_sys_list).T
                df = df.astype(float)


# Read the file using a flexible delimiter (tab or multiple spaces)
#df = pd.read_csv(file_path, na_values='-', index_col=0,header=0, skipinitialspace=True, comment='#')

# Optionally, you can print the dataframe to check the contents
print(df)
df = df*-1

df.index = [r'$Na_2S$', r'$Na_2S_2$', r'$Na_2S_4$', r'$Na_2S_6$', r'$Na_2S_8$']
ax = df.plot(kind='line', figsize=(8, 6), marker='o', linewidth=1.5)

#ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
print(ax.get_xticklabels())

# Rotate individual labels
#for label in ax.get_xticklabels():
#    label.set_rotation(0)


#ax.set_xticklabels([r'$Na_2S$', r'$Na_2S_2$', r'$Na_2S_4$', r'$Na_2S_6$', r'$Na_2S_8$'])
#ax.set_xticklabels(["Al2S3", "Al2S6", "Al2S12", "Al2S18"])
plt.xticks(rotation=0)
plt.axhline(y=0.901966, linestyle='--', linewidth=1, color='red')
# Customize labels if needed
plt.ylabel('Absorption energy of NaPSs (eV)')
plt.ylim(0, 4.5)
#plt.xticks(rotation=90)
#plt.title('Comparison of absorption energy of Na-polysulfides')

# Save as SVG with desired quality
plt.tight_layout()
plt.savefig('absorption_energy_plot_new.svg', format='svg', dpi=300, bbox_inches='tight')

# Optionally display the plot
plt.show()
