#
# Messy but works
# Created by Taylor B. Sage
# pulled from: https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm Javascript and converted into python
# 
# 

import numpy as np

def datum(a, b, f):
    f = 1-b/a
    eccsq = 1 - (b*b)/(a*a)
    ecc = np.sqrt(eccsq)
        
    return {
    "a" : a,
    "b" : b,
    "f" : f,
    "ecc" : ecc,
    "eccsq" : eccsq
    }

#WGS84 Datum
def wgs84Datum():
    wgs84a = 6378.137
    wgs84f = 1.0/298.257223563
    wgs84b = wgs84a * ( 1.0 - wgs84f )
    return datum(wgs84a,wgs84b,wgs84f)

wgs84 = wgs84Datum()


def radcur(lat):

        dtr = np.pi/180.0

        a = wgs84['a']
        b = wgs84['b']

        eccsq  = wgs84['eccsq']
        ecc = wgs84['ecc']

        clat  =  np.cos(dtr*lat)
        slat  =  np.sin(dtr*lat)

        dsq   =  1.0 - eccsq * slat * slat
        d     =  np.sqrt(dsq)

        rn    =  a/d
        rm    =  rn * (1.0 - eccsq ) / dsq

        rho   =  rn * clat
        z     =  (1.0 - eccsq ) * rn * slat
        rsq   =  rho*rho + z*z
        r     =  np.sqrt( rsq )

        return np.array([r, rn, rm])

#convert from Lat lon to cartesian
def latlonelv_xyz(flat, flon, altkmi):
        dtr =  np.pi/180.0

        altkm = altkmi

        clat = np.cos(dtr*flat)
        slat = np.sin(dtr*flat)
        clon = np.cos(dtr*flon)
        slon = np.sin(dtr*flon)

        rrnrm  = radcur (flat)
        rn     = rrnrm[1]
        re     = rrnrm[0]

        ecc    = wgs84['ecc']
        esq    = ecc*ecc

        x =  (rn + altkm) * clat * clon
        y =  (rn + altkm) * clat * slon
        z =  ( (1-esq)*rn + altkm ) * slat

        return  np.array([x, y, z])


def rearth (lati):

    rrnrm =  radcur ( lati )
    r     =  rrnrm[0]
    return r


def gd2gc (flatgd, altkm ):

    dtr   = np.pi/180.0
    rtd   = 1/dtr

    ecc   =  wgs84['ecc']
    esq   =  ecc*ecc

    altnow  =  altkm

    rrnrm   =  radcur (flatgd)
    rn      =  rrnrm[1]

    ratio   = 1 - esq*rn/(rn+altnow)

    tlat    = pi.tan(dtr*flatgd) / ratio
    flatgc  = rtd * pi.arctan(tlat)

    return  flatgc


def gc2gd (flatgc, altkm):

    dtr   = np.pi/180.0
    rtd   = 1/dtr

    ecc   =  wgs84['ecc']
    esq   =  ecc*ecc

    altnow  =  altkm

    rrnrm   =  radcur (flatgc)
    rn      =  rrnrm[1]

    ratio   = 1 - esq*rn/(rn+altnow)

    tlat    = np.tan(dtr*flatgc) / ratio
    flatgd  = rtd * np.arctan(tlat)

    rrnrm   =  radcur ( flatgd )
    rn      =  rrnrm[1]

    ratio   =  1  - esq*rn/(rn+altnow)
    tlat    =  np.tan(dtr*flatgc)/ratio
    flatgd  =  rtd * np.arctan(tlat)

    return  flatgd

#convert from cartesian to lat lon and elv
def xyz_latlonelv ( xvec ):

    # degrees to radians conversion value
    dtr =  np.pi/180.0
    esq    =  wgs84["ecc"]*wgs84["ecc"]

    x      = xvec[0]
    y      = xvec[1]
    z      = xvec[2]

#    rp = np.sqrt(np.dot(xvec,  xvec.T))
    rp = np.sqrt(xvec[0] * xvec[0] +  xvec[1] * xvec[1] +  xvec[2] * xvec[2] )
    
    flatgc = np.arcsin ( xvec[2] / rp )/dtr


    testval = np.abs(xvec[:1]).sum()

    if testval < 1.0e-10:
        flon = 0.0
    else:
        flon = np.arctan2 ( y,x )/dtr

#     if flon < 0.0:
#         flon = flon + 360.0

    p = np.linalg.norm(np.array([x,y]))
    p = np.sqrt(x*x + y*y)
    if p < 1.0e-10:
        flat = 90.0
        if z < 0.0:
            flat = -90.0

        altkm = rp - rearth(flat)
        return  np.array([flat,flon,altkm])

    rnow  =  rearth(flatgc)
    altkm =  rp - rnow
    flat  =  gc2gd (flatgc,altkm)

    rrnrm =  radcur(flat)
    rn    =  rrnrm[1]

    for kount in range(5):
        slat  =  np.sin(dtr*flat)
        tangd =  ( z + rn*esq*slat ) / p
        flatn =  np.arctan(tangd)/dtr

        dlat  =  flatn - flat
        flat  =  flatn
        clat  =  np.cos( dtr*flat )

        rrnrm =  radcur(flat)
        rn    =  rrnrm[1]

        altkm =  (p/clat) - rn

        if np.abs(dlat) < 1.0e-12: 
            break
    
    return  np.array([flat,flon,altkm])

# the elevation of a target (B) above the horizon given in degrees as viewed from pont A
def tgt_elv(a,b):
    b_norm = (b - a) / np.linalg.norm(b - a)
    a_norm = a/np.linalg.norm(a)
    return  90 - np.rad2deg(np.arccos(np.dot(b_norm,a_norm)))


def test_elv():
    # GPS BIIR-11
    # Position:
    # Altitude: 20038 km
    # Lat: 3.1°S
    # Long: 39.2°W

    # From Herndon:
    # Mag: Unknown
    # Dist: 23327 km
    # Alt: 22°
    # Az: 131°

    # 

    #sat position
    target = (-3.1, -39.2,20038)
    #test location position
    test_loc = (38.919380,-77.448290,0)

    target_xyz = latlonelv_xyz(*target)
    test_loc_xyz = latlonelv_xyz(*test_loc)
    print(tgt_elv(test_loc_xyz,target_xyz))

if __name__ == "__main__":
    test_elv()
    x = ''
    x = input()