import utm
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata,interp2d,Rbf
from scipy.fft import fft, ifft
from scipy import stats

#1) DATA COLLECTION

df = pd.read_csv(r"C:\Users\xboxm\vs_workspace\github_projects\gravity_anomaly_splitting\data_grav.csv")

#2) DATA PREPARATION

y = df.LATITUDE
x = df.LONGITUDE
z = df.STATION_ELEV
col = df.BOUGUER_AN

#3 & 4) DATA PROCESSING AND VISUALIZATION

plt.figure(figsize = (12, 10))
plt.subplot(1,2,1)
plt.scatter(x, y, c = z, cmap = "viridis", s=1)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar()
plt.title("Elevation Map (m)")
plt.subplot(1,2,2)
plt.scatter(x, y, c = col, cmap = "jet", s=1)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar()
plt.title("Bouguer Anomaly Map (mGal)")
plt.show()

gridNum = 500
# Min max coordinates
xmin,xmax = min(x),max(x)
ymin,ymax = min(y),max(y)

xgrid = np.linspace(xmin, xmax, gridNum)
ygrid = np.linspace(ymin, ymax, gridNum)
xgrid, ygrid = np.meshgrid(xgrid, ygrid)# Interpolation
cba = griddata((x,y),col,(xgrid,ygrid), method = 'linear')



newxgrid=np.linspace(-4,-1,num=gridNum)
newygrid=np.linspace(58,51,num=gridNum)

newxgrid, newygrid = np.meshgrid(newxgrid, newygrid)

z = griddata((xgrid.reshape(-1),ygrid.reshape(-1)),cba.reshape(-1),(newxgrid,newygrid), method = 'linear')

ip2=z.diagonal()
dx=np.linspace(-4,-1,num=gridNum)
dy=np.linspace(58,51,num=gridNum)

plt.figure(figsize = (6, 10))
plt.scatter(x, y, c = col, cmap = "jet", s=1)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar()
plt.title("Bouguer Anomaly Map with Slice (mGal)")
x1=np.array([-4,-1])
y1=np.array([58,51])
plt.plot(x1, y1)
plt.show()


plt.figure(figsize = (12,5))
plt.subplot(1,2,1)
plt.plot(dx,ip2)
plt.title("Longitude vs CBA (slice)")
plt.xlabel("Longitude")
plt.ylabel("CBA")
plt.subplot(1,2,2)
plt.plot(dy,ip2)
plt.xlim(max(dy),min(dy)) # reverse axis
plt.title("Latitude vs CBA (slice)")
plt.xlabel("Latitude")
plt.ylabel("CBA")
plt.show()


df_slice = pd.DataFrame({"Longitude": dy, "Latitude": dx,"CBA": ip2})

utm_data=utm.from_latlon(dy, dx)
utm_x=np.asarray(utm_data[0])
utm_y=np.asarray(utm_data[1])
df_slice["UTM_X"], df_slice["UTM_Y"] = utm_x,utm_y

# calculate distance
UTM = np.array((df_slice.UTM_X, df_slice.UTM_Y))
xd, yd = np.diff(df_slice.UTM_X), np.diff(df_slice.UTM_Y)
d = np.sqrt(xd**2 + yd**2)
d = np.cumsum(d)
df_slice["dist"] = np.insert(d,0, 0)

plt.figure()
plt.plot(df_slice.dist,ip2)
plt.title("Distance vs CBA")
plt.xlabel("Distance")
plt.ylabel("CBA")
plt.show()

#calculate A
tmp = np.asarray(df_slice.CBA, dtype=float)
A = fft(tmp)
A = np.abs(A)
lnA = np.log(A)
lnA1 = np.log(A**2)

#calculate k
distdiff = df_slice["dist"][len(df_slice)-1]-df_slice["dist"][0]
f = []
fsample=0
for i in range(len(df_slice)):
 fsample = fsample + 1/distdiff
 f.append(fsample)
f = np.array(f)
k = 2*np.pi*f

plt.figure()
plt.plot(k,lnA1)
plt.title("k vs Ln(A^2)")
plt.xlabel('Wavenumber (k)')
plt.ylabel('ln A')
plt.show()


plt.figure()
plt.plot(k[0:len(k)//2],lnA1[0:len(k)//2])
plt.title("k/2 vs Ln(A^2)")
plt.xlabel('Wavenumber (k)')
plt.ylabel('ln A')
plt.show()




reg_cut = 25 #user defined

# Regional cutoff
k_reg_cut = k[0:reg_cut]
lnA_reg_cut = lnA[0:reg_cut]
# Residual cutoff
k_res_cut = k[reg_cut:len(k)//2]
lnA_res_cut = lnA[reg_cut:len(k)//2]

a_reg, b_reg, _, _,std_err1 = stats.linregress(k_reg_cut, lnA_reg_cut)
a_res, b_res,_, _,std_err2= stats.linregress(k_res_cut, lnA_res_cut)
lnA_reg_fit = a_reg*k_reg_cut + b_reg
lnA_res_fit = a_res*k_res_cut + b_res
# Plot spectrum and regressed lines


plt.figure(figsize=(6,5))
plt.scatter(k_reg_cut, lnA_reg_cut, color='b')
plt.plot(k_reg_cut, lnA_reg_fit, '-', color='blue', label='Regional')
plt.scatter(k_res_cut, lnA_res_cut, color='r')
plt.plot(k_res_cut, lnA_res_fit, '-', color='red', label='Residual')
plt.legend()
plt.xlabel('Wavenumber (k)')
plt.ylabel('ln A')
plt.show()


# Calculate cutoff frequency
freq_cut = (b_res - b_reg) / (a_reg - a_res)

# Calculate window
dist = df_slice.dist
window = 2*np.pi/((dist[2]-dist[1])*freq_cut)
window = int(window)

df_slice["Regional"] = df_slice.CBA.rolling(window).mean()
df_slice["Residual"] = df_slice.CBA -df_slice["Regional"]

plt.figure()
plt.plot(df_slice.dist,df_slice.CBA,color='blue', label='Bouguer Anomaly (mGal)')

plt.plot(df_slice.dist,df_slice.Regional,color='black', label='Regional Anomaly (mGal)')

plt.plot(df_slice.dist,df_slice.Residual,color='red', label='Residual Anomaly (mGal)')
plt.title("Anomalies of Slices")
plt.xlabel("Distance (m)")
plt.legend()
plt.show()
