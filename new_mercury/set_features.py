"""
Sets features for line graphs
using mass ratio dictionaries
"""

import pickle as p

mass_ratios = [round(0.05 * x, 2) for x in range(2,11)]

######### VARIABLE SETTINGS BELOW ##########
colors = ['red', 'orange', 'gold', 
          'lime', 'green', 'cyan', 'blue',
          'darkblue', 'purple']

markers = ['D', 'p', 's', '8', 'o', '*', 'v', 'h', '^']
######### VARIABLE SETTINGS ABOVE ##########

# Set up dictionaries
color_dict = {}
marker_dict = {}

for mu,c,m in zip(mass_ratios, colors, markers):
	color_dict[mu] = c
	marker_dict[mu] = m

# Write to files
p.dump(color_dict, open("/Users/Sportsfan77777/planets/plots/features/colors.p", "wb"))
p.dump(marker_dict, open("/Users/Sportsfan77777/planets/plots/features/markers.p", "wb"))

