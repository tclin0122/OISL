# OISL link calculation
import numpy as np
import math
import scipy.constants as constant
# Transmitter parameter
PT_w=0.686*1000#34.05 #mW
PT=10*np.log10(PT_w)
eta_t=0.8 #efficiency
lamb=1550/np.power(10,9) # wavelength m
p_err_t=1/np.power(10,6) #rad
FWHM=15/np.power(10,6) #rad
Gt=16/(np.power(FWHM,2))
GT=10*np.log10(Gt*eta_t)
print(GT)
Lt=np.exp(-Gt*np.power(p_err_t,2))
LT=10*np.log10(Lt)


# Channel parameter
d_ss=4000000 #distance m
L_fsl=20*np.log10(lamb/(4*math.pi*d_ss))
# Performance parameter
data_rate=10*(10**9) #Gbps
BER=1/np.power(10,6) #OCT=10-6

# Receiver parameter
eta_r=0.8

p_err_r=1/np.power(10,6) #rad
D=80/1000 #receiver diameter mm
Gr=np.power((D*math.pi/lamb),2) #receiver gain
GR=10*np.log10(Gr*eta_r)
Lr=np.exp(-Gr*np.power(p_err_r,2))
LR=10*np.log10(Lr)
# should add receiver component loss
k=1.380649*(10**(-23))
T=290 #K
B=data_rate #for OOK bandwidth is data rate
NF=24  #noise figure of receiver <= affect the design of receiver
Rx_s=10*np.log10(k*T*B)+30+11+NF #dBm Thernal noise dBW +30 to dBm +required SNR +NF
Rx_s=-35.5
print('Receiver sensitivity=',Rx_s,'dBm')
#Rx_s=-36.1 #dBm

# Calculation in dB
Pr=PT+GT+GR+LT+LR+L_fsl
Margin=Pr-Rx_s

#BER for OOK-NRZ
snr_linear = 10 ** (Margin / 10.0)
BER=0.5*math.erfc(math.sqrt(snr_linear))


print('EIRP =',PT+GT,'dB')
print('Path Loss =',L_fsl,'dB')
print('Receive power =',Pr,'dB')
print('SNR =',Margin,'dB')
print('BER =',BER)
print('From SDA OCT 3.1.0 OCT-030, Margin should be greater than 3dB')
print("freq=",constant.c/lamb)