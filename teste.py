import math
import numpy
import mpu

#Robson Library
# from PyRadioLoc.Pathloss.Models import FreeSpaceModel
# from PyRadioLoc.Pathloss.Models import FlatEarthModel
# from PyRadioLoc.Pathloss.Models import LeeModel
# from PyRadioLoc.Pathloss.Models import EricssonModel
# from PyRadioLoc.Pathloss.Models import Cost231Model
# from PyRadioLoc.Pathloss.Models import Cost231HataModel
# from PyRadioLoc.Pathloss.Models import OkumuraHataModel
# from PyRadioLoc.Pathloss.Models import Ecc33Model
# from PyRadioLoc.Pathloss.Models import SuiModel

#https://gis.stackexchange.com/questions/66/trilateration-using-3-latitude-longitude-points-and-3-distances/415#415
#https://en.wikipedia.org/wiki/Trilateration
#http://obeattie.github.io/gmaps-radius/?lat=51.500358&lng=-0.125506&z=10&u=mi&r=5

#assuming elevation = 0
earthR = 6371
# LatA = 37.418436
# LonA = -121.963477
# DistA = 0.265710701754
# LatB = 37.417243
# LonB = -121.961889
# DistB = 0.234592423446
# LatC = 37.418692
# LonC = -121.960194
# DistC = 0.0548954278262

# PontoLa = -8.075475
# PontoLo = -34.89629
#Com esses raios absurdos, n ha intersec
#Os circulos tem que overlapear

# LatA = -8.068361111 #BTS1 -8.068361111,-34.892722222 112.99db
# LonA = -34.892722222
# DistA = 0.72 #real eh 0.88km
# DistA = 0.5
#
# LatA = -8.075916667 #BTS2 -8.075916667,-34.894611111 114.39db
# LonA = -34.894611111
# DistA = 1.07 #real eh 0.19km
# DistB = 0.75

# LatB = -8.076361111 #BTS3 -8.076361111,-34.908 128.64db GOOD
# LonB = -34.908
# DistB = 1.72 #real eh 1.29km
# DistC = 1.5

# LatC = -8.075916667 #BTS4 -8.075916667,-34.8946111116 114..64db BTS4 nor BTS1
# LonC = -34.8946111116
# DistC = 0.52086778

# LatB = -8.066 # BTS5 -8.066, -34.8894444444444 128.09db
# LonB = -34.8894444444444
# DistB = 2

# LatC = -8.06458333333333 #BTS6 -8.06458333333333,-34.8945833333333 133.74db GOOD
# LonC = -34.8945833333333
# DistC = 1.78

#using authalic sphere
#if using an ellipsoid this step is slightly different
#Convert geodetic Lat/Long to ECEF xyz
#   1. Convert Lat/Long to radians
#   2. Convert Lat/Long(radians) to ECEF
xA = earthR *(math.cos(math.radians(LatA)) * math.cos(math.radians(LonA)))
yA = earthR *(math.cos(math.radians(LatA)) * math.sin(math.radians(LonA)))
zA = earthR *(math.sin(math.radians(LatA)))

xB = earthR *(math.cos(math.radians(LatB)) * math.cos(math.radians(LonB)))
yB = earthR *(math.cos(math.radians(LatB)) * math.sin(math.radians(LonB)))
zB = earthR *(math.sin(math.radians(LatB)))

xC = earthR *(math.cos(math.radians(LatC)) * math.cos(math.radians(LonC)))
yC = earthR *(math.cos(math.radians(LatC)) * math.sin(math.radians(LonC)))
zC = earthR *(math.sin(math.radians(LatC)))

P1 = numpy.array([xA, yA, zA])
P2 = numpy.array([xB, yB, zB])
P3 = numpy.array([xC, yC, zC])
print P1
print "P2-P1",numpy.linalg.norm(P2-P1)
print "P3-P1",numpy.linalg.norm(P3-P1)
print "P3-P2",numpy.linalg.norm(P3-P2)

#from wikipedia
#transform to get circle 1 at origin
#transform to get circle 2 on x axis
def lateration(P1,P2,P3):
    ex = (P2 - P1)/(numpy.linalg.norm(P2 - P1))
    i = numpy.dot(ex, P3 - P1)


    ey = (P3 - P1 - i*ex)/(numpy.linalg.norm(P3 - P1 - i*ex))
    ez = numpy.cross(ex,ey)
    d = numpy.linalg.norm(P2 - P1)
    j = numpy.dot(ey, P3 - P1)
    print "I", pow(i,2)
    print "J", pow(j,2)
    #from wikipedia
    #plug and chug using above values
    x = (pow(DistA,2) - pow(DistB,2) + pow(d,2))/(2*d)
    y = ((pow(DistA,2) - pow(DistC,2) + pow(i,2) + pow(j,2))/(2*j)) - ((i/j)*x)


    # only one case shown here
    print pow(DistA,2)
    print pow(x,2)
    print pow(y,2)

    #nao pode dar uma raiz negativa
    try:
        z = numpy.sqrt((pow(DistA,2) - pow(x,2) - pow(y,2)))
    except:
        z = float('nan')

    #triPt is an array with ECEF x,y,z of trilateration point
    triPt = P1 + x*ex + y*ey + z*ez

    #convert back to lat/long from ECEF
    #convert to degrees
    lat = math.degrees(math.asin(triPt[2] / earthR))
    lon = math.degrees(math.atan2(triPt[1],triPt[0]))
    #
    print lat, lon

def condition(P1,P2,DistA,DistB):
    # print (numpy.linalg.norm(P2-P1)- DistA)
    # print " < "
    # print DistB
    # print " < "
    # print (numpy.linalg.norm(P2-P1) + DistA)
    # print "////////////"
    if( ((numpy.linalg.norm(P2-P1)- DistA) < DistB)
        and (DistB < (numpy.linalg.norm(P2-P1) + DistA))):
        print "OK"
    else:
        print "NOT OK"

if __name__ == '__main__':
    # condition(P1,P2,DistA,DistB)
    # condition(P1,P3,DistA,DistC)
    # condition(P2,P3,DistB,DistC)
    # lateration(P1,P2,P3)
    # print
    m1 = FreeSpaceModel(1800)
    m2 = FlatEarthModel(1800)
    m3 = LeeModel(1800)
    m4 = EricssonModel(1800)
    m5 =  Cost231Model(1800)
    m6 = Cost231HataModel(1800)
    m7 = OkumuraHataModel(1800)
    m8 = Ecc33Model(1800)
    m9 = SuiModel(2100)
    print("FreeSapce:{}".format(m1.pathloss(0.8)))
    print("FlatEarthModel:{}".format(m2.pathloss([0.8,0.5])))
    print("LeeModel:{}".format(m3.pathloss([0.8,0.5])))
    print("EricssonModel:{}".format(m4.pathloss([0.8,0.5])))
    print("Cost231Model:{}".format(m5.pathloss([0.8,0.5])))
    print("Cost231HataModel:{}".format(m6.pathloss([0.8,0.5])))
    print("OkumuraHataModel:{}".format(m7.pathloss([0.8,0.5])))
    print("Ecc33Model:{}".format(m8.pathloss([0.8,0.5])))
    print("SuiModel:{}".format(m9.pathloss([0.8,0.5])))
#
# condition(P1,P2,DistA,DistB)
# condition(P1,P3,DistA,DistC)
# condition(P3,P2,DistC,DistB)
# condition(P2,P3,DistB,DistC)
