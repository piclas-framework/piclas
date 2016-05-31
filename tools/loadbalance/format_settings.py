from pylab import *
from matplotlib.font_manager import FontProperties
font_size = 9
print ("Font size:", font_size)
font = FontProperties(size = font_size)
rc('font', size=font_size)
# params = {'backend': 'ps',
params = {'axes.labelsize': font_size,
'text.fontsize': font_size,
'legend.fontsize': font_size,
'xtick.labelsize': 8,
'ytick.labelsize': 8}
rcParams.update(params)
#'text.usetex': True,
#'figure.figsize': fig_size}


#rc('font', family="Times New Roman")
#rc('font', family="Roman")
font_family="serif"
rc('font', family=font_family)
print( "Font style:",font_family)
ms = 6


fig_width_pt = 300.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =1.5*fig_width*golden_mean       # height in inches
print ("Figure width:",fig_width,"[in] = ",fig_width*2.54,"[cm])")
print ("Figure height:",fig_height,"[in] = ",fig_height*2.54,"[cm])")
#fig_size = [fig_width,fig_height]

