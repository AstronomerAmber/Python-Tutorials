'''
Plot a Hubble diagram
(python considers tripple-quotes as a super-string that can span several lines
so it's very handy for long comments).
'''

from matplotlib import pyplot as plt
from matplotlib import rcParams   # lets you change defaults
import numpy as np
from astropy.io import ascii

rcParams['font.size'] = 18          # Nice big fonts!
rcParams['font.family'] = 'serif'   # My preference

data = ascii.read('SNIa_DM.dat')

fig,ax = plt.subplots()

# Let's differentiate between CfA (survey==1) and CSP (survey==2) data
labels = ['CfA','CSP']
colors = ['blue','red']
symbols = ['o', 's']

# group the data by survey. I just looked this up in the astropy docs
gdata = data.group_by('survey')

# we can now do a loop to plot the (in this case) 2 surveys
for i in [0,1]:
   ax.errorbar(gdata.groups[i]['zcmb'], gdata.groups[i]['DM'], fmt=symbols[i],
         yerr=gdata.groups[i]['eDM'], capsize=0, mec='black', mfc=colors[i], 
         label=labels[i])
ax.legend(loc='upper left')

# Now let's plot some Hubble laws, one line for each different value of the
# Hubble constant.
zz = np.linspace(data['zcmb'].min(),data['zcmb'].max(), 1000)
for Ho in [50,60,70,80,90,100]:
   # For Ho in km/s/Mpc, the Hubble distance modulus is:
   HDM = 5.0*np.log10(zz*3e5/Ho) + 25
   ax.plot(zz, HDM, '-', color='k')
   # Place a little label at the end of each
   ax.text(zz[-1], HDM[-1], str(Ho), va='center', ha='left', fontsize=10)
ax.set_xlabel('$z_{cmb}$', fontsize=18)
ax.set_ylabel('Distance Modulus', fontsize=16)

# Make space for everything
plt.tight_layout()
fig.savefig('Hubble.pdf')

ax.set_xscale('log')
fig.savefig('Hubble_log.pdf')

