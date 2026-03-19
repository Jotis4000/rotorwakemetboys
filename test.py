import numpy as np

x = np.linspace(0.25,1,50)
xi = [0.25,0.7,1.0]
    # Design Variables
tcoeff = np.polyfit(xi,[1,2,3])
ccoeff = np.polyfit(xi,[5,6,7])

tdist = tcoeff[0]*x**2+tcoeff[1]*x+tcoeff[2]
cdist = ccoeff[0]*x**2+ccoeff[1]*x+ccoeff[2]