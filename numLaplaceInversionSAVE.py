#
#  fabiansThesisCalculations
#  
#
#  Created by Claas Hueter on 5/9/12.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#
from sys import argv, exit
import numpy as np
#import pylab
#from matplotlib.pylab import *

################################################################################
 	
from cmath import *
def Talbot(t):
    N = 24
#   Initiate the stepsize
    h = 2*pi/N;
#   Shift contour to the right in case there is a pole on the positive real axis : Note the contour will
#   not be optimal since it was originally devoloped for function with
#   singularities on the negative real axis
#   For example take F(s) = 1/(s-1), it has a pole at s = 1, the contour needs to be shifted with one
#   unit, i.e shift  = 1. But in the test example no shifting is necessary
    shift = 0.0;
    ans =   0.0;
    
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
        
#   The for loop is evaluating the Laplace inversion at each point theta which is based on the trapezoidal   rule
    for k in range(0,N):
        theta = -pi + (k+1./2)*h;
        z = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans = ans + exp(z*t)*F(z)*dz;
        
    return ((h/(2j*pi))*ans).real        
      

#   Here is the Laplace function to be inverted, should be changed manually        
def F(s):
    return 1.0/(s+1.0)

#   Define the trig functions cot(phi) and csc(phi)
def cot(phi):
    return 1.0/tan(phi)

def csc(phi):
    return 1.0/sin(phi)
#################################################################################       	

def beta(d0_div_delta, phi, Delta):
	return 2*d0_div_delta*phi**2/Delta

def alpha(DeltaPM,Delta):
		if np.abs(Delta) <= 1e-6:
			return 0.
		else:	
			return DeltaPM/Delta

def discriminantOfMuThirdRoot(alphaVal, betaVal ):
	return (0.25*(0.5*np.pi)**(2./3.) * (betaVal-alphaVal-1)**2 + (0.5*np.pi)**(2./3.)*betaVal )**(1./2.)

def muThirdRoot(discriminantVal,alphaVal, betaVal):
	return discriminantVal - 0.5 * (0.5*np.pi)**(1./3.) * (betaVal-alphaVal-1)

def s_star(mu):
	return (2./(np.pi*mu**2))**(1./3.)

def transformOf_dxdy(alphaVal,betaVal,muVal,s):
	#return (muVal-1./s-(alphaVal/(s+betaVal/muVal)))/(muVal*s-(2./np.pi)**(0.5) * s**(-0.5))
	return (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5))
	'''
	if abs(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5)) <= e-6:
		shifted_s = s+0.0001
		X = (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/ (muVal*shifted_s**2-(2./np.pi)**(0.5) * shifted_s**(0.5))  
	else:
		X = (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5))
	return X
	'''
################################################################################	

def Talbot_for_X(mXi,N,alphaVal,betaVal,muVal):
#   Initiate the stepsize
    h = 2*pi/N;
#   Shift contour to the right in case there is a pole on the positive real axis : Note the contour will
#   not be optimal since it was originally devoloped for function with
#   singularities on the negative real axis
#   For example take F(s) = 1/(s-1), it has a pole at s = 1, the contour needs to be shifted with one
#   unit, i.e shift  = 1. But in the test example no shifting is necessary
    shift = muVal+1.
    ans =   0.0;
    
    if mXi < 0:
        print "ERROR:   Inverse transform can not be calculated for -xi < 0"
        return ("Error");
        
#   The for loop is evaluating the Laplace inversion at each point theta which is based on the trapezoidal   rule
    for k in range(0,N):
        theta = -pi + (k+1./2)*h;
        z = shift + N/mXi*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz = N/mXi*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans = ans + exp(z*mXi)*transformOf_dxdy(alphaVal,betaVal,muVal,z)*dz;
        
    return ((h/(2j*pi))*ans).real        
	
################################################################################
 
mXi_array = np.arange(0.01,3.,0.01)
ivantsov_deriv_list = [ (1./(sqrt(2.*float(el)))) for el in list(mXi_array)] 
ivantsov_list = [ sqrt(2.*float(el)) for el in list(mXi_array)]
ivantsov_array = np.asarray(ivantsov_list)
ivantsov_deriv_array = np.asarray(ivantsov_deriv_list)
approx_array = np.asarray([ (1.-2*el/pi + 8*sqrt(2)*el**(1.5)/(3*pi**2)) for el in list(mXi_array)])
#shape_approx_array = np.asarray()


inDelta = 0.05
inDeltaPM_0 = 0.0
inDeltaPM_0p05 = 0.05
inDeltaPM_0p1 = 0.1
inDeltaPM_m0p05 = -0.05
inDeltaPM_m0p1 = -0.1
inPhi = 0.01
ind0_div_delta = 1000


shift_mu = 0.#1.13*e-9. ## trial and error currently for smoothing the singularit at the pole of X(s) - this shifts the singularit away from s_star
s_array = np.arange(0.05,10.,0.05)

##############################################################################
mu1o3Val_0 = muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0,inDelta), beta(ind0_div_delta, inPhi, inDelta))
transformOf_dxdy_array_0 = transformOf_dxdy(alpha(inDeltaPM_0,inDelta),beta(ind0_div_delta, inPhi, inDelta),mu1o3Val_0**3+shift_mu,s_array)
muRes_0 = mu1o3Val_0**3

mu1o3Val_0p05 = muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta))
transformOf_dxdy_array_0p05 = transformOf_dxdy(alpha(inDeltaPM_0p05,inDelta),beta(ind0_div_delta, inPhi, inDelta),mu1o3Val_0p05**3+shift_mu,s_array)
muRes_0p05 = mu1o3Val_0p05**3

mu1o3Val_0p1 = muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta))
transformOf_dxdy_array_0p1 = transformOf_dxdy(alpha(inDeltaPM_0p1,inDelta),beta(ind0_div_delta, inPhi, inDelta),mu1o3Val_0p1**3+shift_mu,s_array)
muRes_0p1 = mu1o3Val_0p1**3

mu1o3Val_m0p05 = muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_m0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_m0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta))
transformOf_dxdy_array_m0p05 = transformOf_dxdy(alpha(inDeltaPM_m0p05,inDelta),beta(ind0_div_delta, inPhi, inDelta),mu1o3Val_m0p05**3+shift_mu,s_array)
muRes_m0p05 = mu1o3Val_m0p05**3

mu1o3Val_m0p1 = muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_m0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_m0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta))
transformOf_dxdy_array_0p1 = transformOf_dxdy(alpha(inDeltaPM_m0p1,inDelta),beta(ind0_div_delta, inPhi, inDelta),mu1o3Val_m0p1**3+shift_mu,s_array)
muRes_m0p1 = mu1o3Val_m0p1**3

print alpha(inDeltaPM_0, inDelta), 'alpha', beta(ind0_div_delta, inPhi, inDelta), 'beta', muRes_0, 'mu'
print alpha(inDeltaPM_0p05, inDelta), 'alpha', beta(ind0_div_delta, inPhi, inDelta), 'beta', muRes_0p05, 'mu'
print alpha(inDeltaPM_0p1, inDelta), 'alpha', beta(ind0_div_delta, inPhi, inDelta), 'beta', muRes_0p1, 'mu'
print alpha(inDeltaPM_m0p05, inDelta), 'alpha', beta(ind0_div_delta, inPhi, inDelta), 'beta', muRes_m0p05, 'mu'
print alpha(inDeltaPM_m0p1, inDelta), 'alpha', beta(ind0_div_delta, inPhi, inDelta), 'beta', muRes_m0p1, 'mu'
##############################################################################
#print transformOf_dxdy_array

filename = 'invTrafo_dxdy_VS_IvantsovSolution_Delta' + str(inDelta) + '_DeltaPM_0' + str(inDeltaPM_0) + '_Phi' + str(inPhi) + '_d_0_div_delta' + str(ind0_div_delta) + '_mu_0_0p05_0p1_m0p05_m0p1' + str(muRes_0) + ' ' + str(muRes_0p05) + ' ' + str(muRes_0p1) + ' ' + str(muRes_m0p05) + ' ' + str(muRes_m0p1) +'.png'
##############################################################################
invTrafoList_0 = [Talbot_for_X( float(el),24,alpha(inDeltaPM_0,inDelta),beta(ind0_div_delta, inPhi, inDelta),muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0,inDelta), beta(ind0_div_delta, inPhi, inDelta))**3+shift_mu) for el in list(mXi_array)]
invTrafo_array_0 = np.asarray( invTrafoList_0 )
mXi_dxdy_List_0 = [ [str(mXi_array[i]), str(invTrafo_array_0[i])]  for i in range(len(invTrafoList_0))]

invTrafoList_0p05 = [Talbot_for_X( float(el),24,alpha(inDeltaPM_0p05,inDelta),beta(ind0_div_delta, inPhi, inDelta),muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta))**3+shift_mu) for el in list(mXi_array)]
invTrafo_array_0p05 = np.asarray( invTrafoList_0p05 )
mXi_dxdy_List_0p05 = [ [str(mXi_array[i]), str(invTrafo_array_0p05[i])]  for i in range(len(invTrafoList_0p05))]

invTrafoList_0p1 = [Talbot_for_X( float(el),24,alpha(inDeltaPM_0p1,inDelta),beta(ind0_div_delta, inPhi, inDelta),muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta))**3+shift_mu) for el in list(mXi_array)]
invTrafo_array_0p1 = np.asarray( invTrafoList_0p1 )
mXi_dxdy_List_0p1 = [ [str(mXi_array[i]), str(invTrafo_array_0p1[i])]  for i in range(len(invTrafoList_0p1))]

invTrafoList_m0p05 = [Talbot_for_X( float(el),24,alpha(inDeltaPM_m0p05,inDelta),beta(ind0_div_delta, inPhi, inDelta),muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_m0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_m0p05,inDelta), beta(ind0_div_delta, inPhi, inDelta))**3+shift_mu) for el in list(mXi_array)]
invTrafo_array_m0p05 = np.asarray( invTrafoList_m0p05 )
mXi_dxdy_List_m0p05 = [ [str(mXi_array[i]), str(invTrafo_array_m0p05[i])]  for i in range(len(invTrafoList_m0p05))]

invTrafoList_m0p1 = [Talbot_for_X( float(el),24,alpha(inDeltaPM_m0p1,inDelta),beta(ind0_div_delta, inPhi, inDelta),muThirdRoot(discriminantOfMuThirdRoot(alpha(inDeltaPM_m0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta) ),alpha(inDeltaPM_m0p1,inDelta), beta(ind0_div_delta, inPhi, inDelta))**3+shift_mu) for el in list(mXi_array)]
invTrafo_array_m0p1 = np.asarray( invTrafoList_m0p1 )
mXi_dxdy_List_m0p1 = [ [str(mXi_array[i]), str(invTrafo_array_m0p1[i])]  for i in range(len(invTrafoList_m0p1))]
##############################################################################


shape_List_0 = []
shape_List_0p05 = []
shape_List_0p1 = []
shape_List_m0p05 = []
shape_List_m0p1 = []

for i in range(len(mXi_array)):
	var_Xi_array = mXi_array[0:i]
	var_invTrafo_array_0 = invTrafo_array_0[0:i]
	var_invTrafo_array_0p05 = invTrafo_array_0p05[0:i]
	var_invTrafo_array_0p1 = invTrafo_array_0p1[0:i]
	var_invTrafo_array_m0p05 = invTrafo_array_m0p05[0:i]
	var_invTrafo_array_m0p1 = invTrafo_array_m0p1[0:i]
	shape_List_0.append(np.trapz(var_invTrafo_array_0, var_Xi_array))
	shape_List_0p05.append(np.trapz(var_invTrafo_array_0p05, var_Xi_array))
	shape_List_0p1.append(np.trapz(var_invTrafo_array_0p1, var_Xi_array))
	shape_List_m0p05.append(np.trapz(var_invTrafo_array_m0p05, var_Xi_array))
	shape_List_m0p1.append(np.trapz(var_invTrafo_array_m0p1, var_Xi_array))

#print shape_List

shape_array_0 = np.asarray(shape_List_0)
shape_array_0p05 = np.asarray(shape_List_0p05)
shape_array_0p1 = np.asarray(shape_List_0p1)
shape_array_m0p05 = np.asarray(shape_List_m0p05)
shape_array_m0p1 = np.asarray(shape_List_m0p1)
#shapeAt_xi = np.trapz(invTrafo_array, mXi_array)
#print shapeAt_xi

#print mXi_dxdy_List

filename_dat = 'resultShapes5Alphas.dat'
#filename_dat = 'invTrafo_dxdy_VS_IvantsovSolution_Delta' + str(inDelta) + '_DeltaPM' + str(inDeltaPM) + '_Phi' + str(inPhi) + '_d_0_div_delta' + str(ind0_div_delta) + '_mu' + str(muRes) + '.dat'
print filename_dat, "is filename_dat"

fdat = open(filename_dat, 'w')
for el in mXi_dxdy_List_0:
	fdat.write((str(el).split())[0].strip('[').strip("'").rstrip(",").rstrip("'") + ' ' + (str(el).split())[1].strip("'").strip(']').rstrip("'") +'\n')
fdat.close()

#mXi_array, invTrafo_array_0,mXi_array, invTrafo_array_0p05,mXi_array, invTrafo_array_0p1,mXi_array,ivantsov_deriv_array, mXi_array, approx_array,
#fig = figure()
#rc('text', usetex=True)
#xlabel(r'-$\xi$')
#ylabel(r'x($\xi$)')
#plot( mXi_array, shape_array_0, mXi_array, shape_array_0p05, mXi_array, shape_array_0p1, mXi_array, shape_array_m0p05, mXi_array, shape_array_m0p1, mXi_array, ivantsov_array)
#title('')
#grid(True)
#fig.savefig('result5Shapes.pdf')

#pylab.show()

'''
from matplotlib.pylab import *
fig = figure()
# ...
# Various commands to generate the plot
# ...
fig.savefig('name_of_plot.png')
'''

'''
fig.savefig('name_of_plot.pdf',
             dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait',
             papertype='letter')
'''

'''
pylab.plot(data) # pylab will automatically add an x-axis from 0 to len(data) - 1

# first argument is when ticks should appear on the x-axis
# the second argument is what the label for each tick should be
# same as -> pylab.xticks([0, 60, 120...], [0, 1, 2...])
pylab.xticks(range(0, len(data), 60), range(len(data)/60))

# let the reader know the units
pylab.xlabel("hours")
'''
'''
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
# plot figure
# ...
# annotate figure
plt.xlabel(r'$\mu$ = 50')
plt.ylabel(r'$\sigma$ = 1.5')
'''
'''
import pylab
pylab.rc('text', usetex=True)
pylab.xlabel(r'my data $\alpha$')
pylab.plot(range(0,5))
pylab.show()
'''


















	
	
