import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
import scipy.ndimage as ndimage


#Compute a rupture Velocity map from the contatenated RF file
#the fault has to be parametrizable in function of x and z
#in fact t and y are mapped in function of the structured grid xi, zi 
#in order to apply the gradient function of numpy


#####Parameters###########

cplotfile="ParRF25-100-RF-concat.dat"
cplotfile="ParRF10-25-RF-concat.dat"
cplotfile="tpv29-50p-RF-concat.dat"
#cplotfile="tpv29Z-all-RF-concat.dat"

#cplotfile="ParRF100-600-RF-concat.dat"
#cplotfile="ParRF250-RF-concat.dat"
cplotfile="/export/data/ulrich/NewRoughFault/RFfiles/NRF10-400-RF-concat-small.dat"

GDmethod='cubic'
GDmethod='linear'
#GDmethod='nearest'
# Size of regular grid
dx, dz = 50., 50.
#use Gaussian Filter?
useGaussianFilter = True
#min max percentiles of V for plotting (eliminate artefacts)
p1, p2=1, 99
#plot Vr=f(x) and Vr=f(y)
plotHist = False
#boundaries values for Vr plot
xm,xp = -20e3, 20e3
zm,zp = -20e3, 0
##########################

#Read Cplot file
#xyt = np.loadtxt(cplotfile, skiprows=20)
xyt = np.loadtxt(cplotfile,  skiprows=1)
print('done reading %s' %cplotfile)
x = xyt[:,0]
y = xyt[:,1]
z = xyt[:,2]
t = xyt[:,3]


# Generate a regular grid to interpolate the data.
Xi = np.arange(min(x), max(x), dx)
Zi = np.arange(min(z), max(z), dz)
xi, zi = np.meshgrid(Xi, Zi)
print('done meshgrid')

# Interpolate using delaunay triangularization 
ti = griddata((x,z), t, (xi, zi), method=GDmethod, fill_value = 1e20)
yi = griddata((x,z), y, (xi, zi), method=GDmethod, fill_value = 0.)
print('done griddata')


if useGaussianFilter:
   # Increase the value of sigma to increase the amount of blurring.
   # order=0 means gaussian kernel
   ti = ndimage.gaussian_filter(ti, sigma=2.0, order=0)
   yi = ndimage.gaussian_filter(yi, sigma=2.0, order=0)

grad = np.gradient(ti)
gradY = np.gradient(yi, dx, dz)

print('done gradient')

dy1 = gradY[0]
dy2 = gradY[1]

# 2d numpy array are indexed as row, column meaning that we should write dy, dx = np.gradient(zi)
# and not dx, dy = np.gradient(zi). That why the indexes are inverted
dtdx=grad[1]/np.sqrt(pow(dx,2)+pow(dy2,2))
dtdz=grad[0]/np.sqrt(pow(dx,2)+pow(dy1,2))

slowness = np.sqrt(np.square(dtdx) + np.square(dtdz))
V=1./ slowness
print('done Vr')

#process data where NaN
#where_are_NaNs = np.isnan(V)
#V[where_are_NaNs] = 0
#where_are_null = np.where(V<1.)
#V[where_are_null] = np.nan

#Show a few percentiles for helping setting up p1 and p2
gradV = np.gradient(V)
dVdx=gradV[1]/np.sqrt(pow(dx,2)+pow(dy2,2))
dVdz=gradV[0]/np.sqrt(pow(dx,2)+pow(dy1,2))

normGV = np.sqrt(np.square(dVdx)+np.square(dVdz))

if False:
	for i in range(1,20):
	   V1=np.percentile(V, i)
	   print("percentile %d: %f" %(i,V1))

	print(" ")
	for i in range(99,70,-1):
	   V1=np.percentile(V, i)
	   print("percentile %d: %f" %(i,V1))

V1=np.percentile(V, p1)
V50=np.percentile(V, 50)
V2=np.percentile(V, p2)

#print("percentiles %d, 50, %d: %d %d %d" %(p1,p2,V1,V50,V2))

#plt.hist(V, bins=[1000*i for i in range(0,10)])
#plt.show()

# Plot the results

fig = plt.figure()
fig.add_subplot(211)

masked_array=np.ma.masked_where(V>1e10, V)
cmap = cm.jet
cmap.set_bad('w',1.)
plt.pcolormesh(xi,zi,V)
#Eliminate Vr artefacts (at rupture rim)
plt.pcolormesh(xi,zi,masked_array,cmap=cmap)
#plt.clim(V1,V2)
plt.clim(0,5400)
plt.colorbar()


CS = plt.contour(Xi, Zi, ti,range(1,21),colors='k')
plt.clabel(CS, fontsize=9, inline=1, fmt='%d')
plt.xlim(xm,xp)
plt.ylim(zm,zp)
plt.axis('equal')
plt.grid()

plotGradientField = False
if plotGradientField:
   downsampling=30
   xi2 = xi[::downsampling,::downsampling]
   zi2 = zi[::downsampling,::downsampling]
   if False:
      grad2a = grad[0][::downsampling,::downsampling]
      grad2b = grad[1][::downsampling,::downsampling]
      normgrad = np.sqrt(grad2a**2 + grad2b**2)
      grad2an = grad2a/normgrad
      grad2bn = grad2b/normgrad
      #suppress some artifacts
      grad2an[0,:] = 0
      grad2bn[0,:] = 0
      grad2an[:,0] = 0
      grad2bn[:,0] = 0
      plt.quiver(xi2, zi2, grad2bn,grad2an,
		 color='0.5', scale=60, width=0.002)

   grad2a = gradY[0][::downsampling,::downsampling]
   grad2b = gradY[1][::downsampling,::downsampling]
   normgrad = np.sqrt(grad2a**2 + grad2b**2)
   normgrad30 = np.percentile(normgrad,30)
   grad2an = grad2a/normgrad30
   grad2bn = grad2b/normgrad30
   normgrad30 = np.percentile(grad2an,30)
   print normgrad30
   normgrad30 = np.percentile(grad2bn,30)
   print normgrad30
   #suppress some artifacts
   grad2an[0,:] = 0
   grad2bn[0,:] = 0
   grad2an[:,0] = 0
   grad2bn[:,0] = 0
   plt.quiver(xi2, zi2, grad2bn,grad2an,
              color='k', scale=100, width=0.002)




#plt.show()
fig.add_subplot(212)

# Plot the results

#normgrad = np.sqrt(dtdx**2 + dtdz**2)
#SlopeInGradDir = (np.multiply(dtdx,dy2)+np.multiply(dtdz,dy1))/normgrad
normgrad = np.sqrt(grad[0]**2 + grad[1]**2)
indexes = np.where(normgrad==0)
normgrad[indexes]=1.0
#SlopeInGradDir = (np.multiply(grad[1],dy2)+np.multiply(grad[0],dy1))/normgrad
#Magnitude of gradient rather than SlopeInGradDir
SlopeInGradDir = np.sqrt(gradY[0]**2 + gradY[1]**2)
#SlopeInGradDir = np.norm(gradY)

print np.percentile(SlopeInGradDir, 20)
print np.percentile(SlopeInGradDir, 50)
print np.percentile(SlopeInGradDir, 80)
#normGradient = np.sqrt(np.square(dy1) + np.square(dy2))

#normR = np.sqrt(np.square(xi-5000)+np.square(zi+10000))
#invnormR = 1/normR
#gradientRadial = ((xi-5000)*dy1+(zi+10000)*dy2) * invnormR

#V222= np.sqrt(np.square(1/dtdx) + np.square(1/dtdz))


plt.pcolormesh(xi,zi,SlopeInGradDir, cmap=cmap)
plt.clim(0.0,0.6)
#plt.pcolormesh(xi,zi-2e4, gradientRadial, cmap=cmap)
#plt.pcolormesh(xi,zi-2e4, V222, cmap=cmap)
#plt.clim(V1,V2)
plt.colorbar()
CS = plt.contour(Xi, Zi, ti,range(1,21),colors='k')
plt.clabel(CS, fontsize=9, inline=1, fmt='%d')

#CS2 = plt.contour(Xi, Zi, normGV, [10], colors='r', linewidths=[2])
#CS3 = plt.contour(Xi, Zi-2e4, V, [1500], colors='b')
#plt.clabel(CS2, fontsize=9, inline=1, fmt='%d')




#plt.xlim(xm,xp)
plt.ylim(zm,zp)
plt.axis('equal')
#plt.clim(0,100)

#ax.set_yticks(numpy.arange(-2e4,0.,2e3))
plt.grid()
plt.show()

plt.figure()

#plt.plot(normGradient.flatten(),V.flatten(),'k.',alpha=0.1)
#plt.plot(gradientRadial.flatten(),V.flatten(),'k.',alpha=0.1)
#plt.ylim(V1,V2)
#plt.show()



if plotHist:
	n=np.shape(V)[0]
	lV1, lV50, lV2 = np.zeros(n), np.zeros(n), np.zeros(n)
	for i in range(n):
	    subV=V[i,:]
	    lV1[i]=np.percentile(subV, 33)
	    lV50[i]=np.percentile(subV, 50)
	    lV2[i]=np.percentile(subV, 67)
	plt.plot(lV1, Zi, label = '33%')
	plt.plot(lV50, Zi, label = '50%')
	plt.plot(lV2, Zi, label = '67%')
	plt.legend()
	plt.title("Vr=f(y)")
	plt.show()

	n=np.shape(V)[1]
	lV1, lV50, lV2 = np.zeros(n), np.zeros(n), np.zeros(n)
	for i in range(n):
	    subV=V[:,i]
	    lV1[i]=np.percentile(subV, 33)
	    lV50[i]=np.percentile(subV, 50)
	    lV2[i]=np.percentile(subV, 67)
	plt.plot(Xi, lV1, label = '33%')
	plt.plot(Xi, lV50, label = '50%')
	plt.plot(Xi, lV2, label = '67%')
	plt.legend()
	plt.title("Vr=f(x)")
	plt.show()

