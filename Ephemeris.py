#initialize orekit and JVM
import orekit_jpype as orekit
orekit.initVM()
# setup the orekit data loading, the file orekit-data.zip shall be in same directory as notebook.
from orekit_jpype.pyhelpers import setup_orekit_curdir, download_orekit_data_curdir
from jpype import JImplements, JOverride
# download_orekit_data_curdir()
setup_orekit_curdir()
# import all the library that used
from org.orekit.time import TimeScalesFactory
from org.orekit.time import AbsoluteDate
from org.orekit.utils import Constants
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle import TLEPropagator
from org.hipparchus.geometry.euclidean.threed import Line
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
from org.orekit.propagation.sampling import OrekitFixedStepHandler
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateComponents, TimeComponents
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import math
from IPython.display import FileLink, display
import nums_from_string
import scipy.constants as constants

# start calulation
#set up time system
utc = TimeScalesFactory.getUTC()


# import the Two Line Elements
# STARLINK-2232
TLE_LINE1="1 48578U 21041AB  23194.63168809 -.00000388  00000+0 -71564-5 0  9998"
TLE_LINE2="2 48578  53.0546 166.8878 0001930  85.3045 274.8164 15.06390900119432"
# STARLINK-2191 
tle_LINE1="1 48565U 21041N   23194.63066506 -.00000644  00000+0 -24326-4 0  9993"
tle_LINE2="2 48565  53.0545 166.3746 0001744  93.1087 267.0101 15.06400969119445"


# Simulation parameters
DURATION = 1.0 * 60 * 60 * 1.5
STEP = 10.0
ANGULAR_OFFSET = 35 # Sensor half width

# setup the propagator
TimeScalesFactory

tle1 = TLE(TLE_LINE1, TLE_LINE2)
tle2 = TLE(tle_LINE1, tle_LINE2)

propagator1 = TLEPropagator.selectExtrapolator(tle1)
propagator2 = TLEPropagator.selectExtrapolator(tle2)
start = tle1.getDate() 

#Create a custom step handler
@JImplements(OrekitFixedStepHandler)
class CorridorHandler():
    
    def __init__(self, angle):
            # set up Earth model
            self.earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                         Constants.WGS84_EARTH_FLATTENING,
                                         FramesFactory.getITRF(IERSConventions.IERS_2010, False))
            
            # set up position offsets, using Earth radius as an arbitrary distance
            self.deltaR = Constants.WGS84_EARTH_EQUATORIAL_RADIUS * math.cos(math.radians(angle))
            self.deltaC = Constants.WGS84_EARTH_EQUATORIAL_RADIUS * math.sin(math.radians(angle))

            # prepare an empty corridor
            self.dates = []
            self.lefts = []
            self.centers = []
            self.rights = []
    
    @JOverride  
    def init(self, s0, t, step):
        # needs to be stated to fulfill the interface specification
        pass

    @JOverride  
    def handleStep(self, currentState):
        # compute sub-satellite track
        date    = currentState.getDate()
        pvInert = currentState.getPVCoordinates()
        t       = currentState.getFrame().getTransformTo(self.earth.getBodyFrame(), date)
        p       = t.transformPosition(pvInert.getPosition())
        v       = t.transformVector(pvInert.getVelocity())
        center  = self.earth.transform(p, self.earth.getBodyFrame(), date)
        
        # compute left and right corridor points
        nadir      = p.normalize().negate()
        crossTrack = p.crossProduct(v).normalize()
        leftLine   = Line(p, Vector3D(1.0, p, self.deltaR, nadir,  self.deltaC, crossTrack), 1.0)
        left       = self.earth.getIntersectionPoint(leftLine, p, self.earth.getBodyFrame(), date)
        rightLine  = Line(p, Vector3D(1.0, p, self.deltaR, nadir, -self.deltaC, crossTrack), 1.0)
        right      = self.earth.getIntersectionPoint(rightLine, p, self.earth.getBodyFrame(), date)
        

        # add the corridor points
        self.dates.append(date)
        self.lefts.append(left)
        self.centers.append(center)
        self.rights.append(right)
        
    @JOverride
    def finish(self, s):
        pass
handler = CorridorHandler(ANGULAR_OFFSET)
handler2 = CorridorHandler(ANGULAR_OFFSET)
propagator1.setStepHandler(STEP, handler)
propagator1.propagate(start, start.shiftedBy(DURATION));
propagator2.setStepHandler(STEP, handler2)
propagator2.propagate(start, start.shiftedBy(DURATION));

def geoline(geopoints):
    lon = [math.degrees(x.getLongitude()) for x in geopoints]
    lat = [math.degrees(x.getLatitude()) for x in geopoints]
    return lon, lat


ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()

lon, lat = geoline(handler.centers)
ax.plot(lon, lat, transform=ccrs.Geodetic(), alpha=0.6, color='green', zorder=3);

lon, lat = geoline(handler2.centers)
ax.plot(lon, lat, transform=ccrs.Geodetic(), alpha=0.6, color='blue', zorder=3);


initialDate=tle1.getDate()
pvs1 = []
pvs2 = []
extrapDate = initialDate
finalDate = extrapDate.shiftedBy(60.0 * 60 * 1 *1)  #seconds

inertialFrame = FramesFactory.getEME2000()

while (extrapDate.compareTo(finalDate) <= 0.0):
    pv1 = propagator1.getPVCoordinates(extrapDate, inertialFrame)
    pvs1.append(pv1)
    pv2 = propagator2.getPVCoordinates(extrapDate, inertialFrame)
    pvs2.append(pv2)
    extrapDate = extrapDate.shiftedBy(20.0)

prop_data1 = pd.DataFrame(data=pvs1, columns=['pv'])

prop_data1['Position'] = prop_data1['pv'].apply(lambda x: x.getPosition())
prop_data2 = pd.DataFrame(data=pvs2, columns=['pv'])
prop_data2['Position'] = prop_data2['pv'].apply(lambda x: x.getPosition())

# fetch data from the propagation data
def fetch_data(pvs1,pvs2,timetag):
    adata=nums_from_string.get_nums(str(pvs1[timetag]))
    bdata=nums_from_string.get_nums(str(pvs2[timetag]))
    (X1,Y1,Z1)=(adata[6],adata[7],adata[8])
    (X2,Y2,Z2)=(bdata[6],bdata[7],bdata[8])
    (Vx1,Vy1,Vz1)=(adata[9],adata[10],adata[11])
    (Vx2,Vy2,Vz2)=(bdata[9],bdata[10],bdata[11])
    print(X1,Y1,Z1,X2,Y2,Z2)
    print(Vx1,Vy1,Vz1,Vx2,Vy2,Vz2)
    return X1,X2,Y1,Y2,Z1,Z2,Vx1,Vy1,Vz1,Vx2,Vy2,Vz2
#calculate direction angle
def dir(X1,Y1,Z1,X2,Y2,Z2):
    sat_1=np.array([[X1], [Y1], [Z1]])
    sat_2=np.array([[X2], [Y2], [Z2]])
    s=sat_2-sat_1
    dist=np.sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2])
    if(s[1]==0 and s[0]>0): azimuth=0
    elif(s[1]==0 and s[0]<0): azimuth=180
    elif(s[0]==0 and s[1]>0): azimuth=90
    elif(s[0]==0 and s[1]<0): azimuth=270
    elif(s[1]>0 and s[0]>0): azimuth=np.rad2deg(math.atan(s[0]/s[1]))
    elif(s[1]>0 and s[0]<0):azimuth=np.rad2deg(math.atan(s[0]/s[1]))+180
    elif(s[1]<0 and s[0]<0):azimuth=np.rad2deg(math.atan(s[0]/s[1]))+180
    elif(s[1]<0 and s[0]>0):azimuth=np.rad2deg(math.atan(s[0]/s[1]))+360
    else: azimuth=0
    if(dist==0):
        print('Satellite collide')
        elevation=0
    else:
        zenith=np.rad2deg(math.acos(s[2]/dist))
        elevation=90-zenith
    return azimuth, elevation, dist
def pointing_ahead(Vx1,Vy1,Vz1,Vx2,Vy2,Vz2):
    dV=np.sqrt(np.power((Vx2-Vx1),2)+np.power((Vy2-Vy1),2)+np.power((Vz2-Vz1),2)) #km
    ahead_angle=2*dV/(constants.c)
    print(dV,'m/s')
    print(ahead_angle*(10**6),"urad")
    return dV, ahead_angle
def doppler_effect(dV,lamb):
    df=dV/lamb
    print('Doppler shift=',df/(10**(9)),'GHz')
    return df

timetag=180
print(len(prop_data1))
(X1,Y1,Z1,X2,Y2,Z2,Vx1,Vy1,Vz1,Vx2,Vy2,Vz2)=fetch_data(pvs1,pvs2,timetag)
az,el,s=dir(X1,Y1,Z1,X2,Y2,Z2)
dV,ahead_angle=pointing_ahead(Vx1,Vy1,Vz1,Vx2,Vy2,Vz2)
print(az,el,s)
lamb=1550*(10**(-9))
df=doppler_effect(dV,lamb)