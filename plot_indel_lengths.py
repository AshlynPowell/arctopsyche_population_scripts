import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullLocator
import numpy as np

max_len = 100
# max_len = 1213 # mafft
# max_len = 1188 # muscle
# max_len = 2201 # prank

lengths = []

# with open("indels_all.csv") as file:
with open("indels_all_muscle.csv") as file:
# with open("indels_all_prank.csv") as file:
    file.readline()
    for line in file:
        items = line.strip().split(",")
        if len(items) == 1:
            break
        
        lengths.append(int(items[5]))
        
all_counts = []
for i in range(1,max_len+1):
    all_counts.append(lengths.count(i))

width = 0.5

fig, ax = plt.subplots()
fig.set_size_inches(8,3)
bottom = np.zeros(max_len)
ax.bar(list(range(1,max_len+1)), all_counts, width, bottom=bottom, color="#020080")
ax.set_xlabel("Indel Length (Amino Acids)", size=12)
ax.set_ylabel("Count", size=12)
ax.tick_params(axis='both', which='major', labelsize=10)
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(NullLocator())
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(25))
ax.margins(x=0.01)
ax.spines[['right', 'top']].set_visible(False)
# plt.savefig("indel_lengths(100).pdf",format='pdf',dpi=1200,bbox_inches='tight', pad_inches=0.25)
plt.savefig("indel_lengths(100)_muslce.pdf",format='pdf',dpi=1200,bbox_inches='tight', pad_inches=0.25)
# plt.savefig("indel_lengths(100)_prank.pdf",format='pdf',dpi=1200,bbox_inches='tight', pad_inches=0.25)


