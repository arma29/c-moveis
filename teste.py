import math
import numpy
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

LatA = -8.068361111 #BTS1 -8.068361111,-34.892722222 112.99db
LonA = -34.892722222
# DistA = 5.91190 #real eh 0.88km
DistA = 0.6
#
LatB = -8.075916667 #BTS2 -8.075916667,-34.894611111 114.39db
LonB = -34.894611111
# DistB = 6.94588 #real eh 0.19km
DistB = 0.7
#
LatC = -8.076361111 #BTS3 -8.076361111,-34.908 128.64db
LonC = -34.908
# DistC = 35.82840 #real eh 1.29km
DistC = 3.5

# LatC = -8.075916667 #BTS4
# LonC = -34.8946111116
# DistC = 7.06687

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
    z = numpy.sqrt((pow(DistA,2) - pow(x,2) - pow(y,2)))

    #triPt is an array with ECEF x,y,z of trilateration point
    triPt = P1 + x*ex + y*ey + z*ez

    #convert back to lat/long from ECEF
    #convert to degrees
    lat = math.degrees(math.asin(triPt[2] / earthR))
    lon = math.degrees(math.atan2(triPt[1],triPt[0]))
    #
    print lat, lon

def condition(P1,P2,DistA,DistB):
    print (numpy.linalg.norm(P2-P1)- DistA)
    print " < "
    print DistB
    print " < "
    print (numpy.linalg.norm(P2-P1) + DistA)
    print "////////////"
    if( ((numpy.linalg.norm(P2-P1)- DistA) < DistB)
        and (DistB < (numpy.linalg.norm(P2-P1) + DistA))):
        print "OK"
    else:
        print "NOT OK"

if __name__ == '__main__':
    condition(P1,P2,DistA,DistB)
    lateration(P1,P2,P3)
#
# condition(P1,P2,DistA,DistB)
# condition(P1,P3,DistA,DistC)
# condition(P3,P2,DistC,DistB)
# condition(P2,P3,DistB,DistC)
