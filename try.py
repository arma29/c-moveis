#!/usr/bin/python
#Similar to TOA
#Try to write a program to solve the exercise of Trilateration GO
#useful http://www.meridianoutpost.com/resources/etools/calculators/calculator-latitude-longitude-distance.php?
import csv
import numpy as np
import math
import mpu
# import pyproj as pj


#Constants
data_BTS = csv.reader(open('Dados/dados_BTSs.csv', 'rb'))
data_LOC = csv.reader(open('Dados/LocTreino_Equipe_4.csv', 'rb'))
# 1501 linhas -> 1500 pontos

lat_index = 1
lon_index = 2
d0 = 0.1	#must be greater than 0
PL = 2 #Pl expoent in free space
B = 6 #number of BTSs
L = 3 #Dimension 2D
P = 1500 #number of points
GSM = 1800 #frequency for FPSL MHZ and meters (-27.55)https://en.wikipedia.org/wiki/Free-space_path_loss

# LOC_mtz = np.zeros([P,L]) #Real points
LOC_mtz = np.array([[1.]*L]*P)
# BTS_mtz = np.zeros([B,L]) #Bts Points
BTS_mtz = np.array([[1.]*L]*B)

#TODO: Decisions and completness test
# def gps_ecef(lat, lon, alt):
# 	ecef = pj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
# 	lla = pj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
# 	x, y, z = pj.transform(lla, ecef, lon, lat, alt, radians=False)

# 	return np.array([x,y,z])

def coord_conv(lat,lon,alt):
	lat = float(lat)
	lon = float(lon)
	earthR = 6371
	x = earthR *(math.cos(math.radians(lat)) * math.cos(math.radians(lon)))
	y = earthR *(math.cos(math.radians(lat)) * math.sin(math.radians(lon)))
	z = earthR *(math.sin(math.radians(lat)))

	return np.array([x,y,z])

def convert_np():
	list_bts = list(data_BTS)
	list_bts.pop(0)

	for i in range(0,B):
		# BTS_mtz[i] = gps_ecef(list_bts[i][1], list_bts[i][2], 0)
		BTS_mtz[i] = coord_conv(list_bts[i][1], list_bts[i][2], 0)
		# for j in range (0,L):
		# 	BTS_mtz[i][j] = list_bts[i][j+1]

	print BTS_mtz


def ML_estimate(vi):
	return d0*pow(10 , (float(vi)/(10*PL)))

def FPSL_estimate(vi):
	#FSPL(db), d in meters, f in MHZ , cte = -27.55
	#FSPL(db), d in meters, f in GHZ , cte = -87.55
	#FSPL(db), d in Kmeters, f in MHZ , cte = +32.45
	#FSPL(db), d in Kmeters, f in GHZ , cte = +92.45
	factor = 1
	f = 1800 #original MHZ
	cte = +32.55
	return math.sqrt(math.sqrt(pow(10, (( float(vi)-cte)-(20*math.log(f,10)))/20)))

#return distance to target array and modify the real points
def dist_array():

	#Read .csv file with 1500 points
	loc_list = list(data_LOC)
	loc_list.pop(0)

	dist_mtz = np.array([[1.]*B]*P)

	#Each point will a array of 6 distance (anchor - target)
	for i in range(0,P):
		# LOC_mtz[i] = gps_ecef(loc_list[i][1], loc_list[i][2],0) #lat, lon ,hei
		LOC_mtz[i] = coord_conv(loc_list[i][1], loc_list[i][2],0)
		for j in range (0,L):
			LOC_mtz[i][j] = loc_list[i][j+1]
		for j in range (0,B):
			# print loc_list[i][j+3]
			dist_mtz[i][j]  = ML_estimate((loc_list[i][j+3]))

			#pathBTS1 to pathBTS6
			# dist_mtz[i][j] = FPSL_estimate(loc_list[i][j+3])


	# print "Dist ponto 1 -",  dist_mtz[0]
	# print "Real ponto 1 -", LOC_mtz[0]
	return dist_mtz

def coefficient_A(i,j,k,m):
	A = np.array([[1.]*L]*L)
	# aux = np.zeros([L,L])
	aux = np.array([BTS_mtz[i],]*L) #mtz de x1*L

	A[0] = np.copy(BTS_mtz[j])
	A[1] = np.copy(BTS_mtz[k])
	A[2] = np.copy(BTS_mtz[m])

	# aux[0] = np.copy(BTS_mtz[i])
	# aux[1] = np.copy(BTS_mtz[i])
	# aux[2] = np.copy(BTS_mtz[i])

	A = np.subtract(A,aux)

	return A

def coefficient_b(i,j,k,m):
	b = np.zeros([L])

	xi = pow(np.linalg.norm(BTS_mtz[i]),2)
	xj = pow(np.linalg.norm(BTS_mtz[j]),2)
	xk = pow(np.linalg.norm(BTS_mtz[k]),2)
	xm = pow(np.linalg.norm(BTS_mtz[m]),2)

	b[0] = xj - xi
	b[1] = xk - xi
	b[2] = xm - xi

	return b

#TODO: new mtz[(len(loc_list),L)] , contendo todos os pontos predictos
def rss_lat(i_index,j_index,k_index,m_index):
	# Ax = b -> x = Ainv*b
	# BTS_mtz = convert_np() #6 pontos das BTS, pode ser fixo tambem

	convert_np() #initializate BTS_mtz
	A = coefficient_A(i_index,j_index,k_index, m_index)

	dist_mtz = dist_array()

	x = np.zeros([P,L])
	for i in range (0,1): #3 to P
		b = coefficient_b(i_index,j_index,k_index, m_index)

		di = pow(dist_mtz[i][i_index],2)
		dj = pow(dist_mtz[i][j_index],2)
		dk = pow(dist_mtz[i][k_index],2)
		dm = pow(dist_mtz[i][m_index],2)

		b[0] = b[0] - (dj - di)
		b[1] = b[1] - (dk - di)
		b[2] = b[2] - (dm - di)

		b = b*1/2
		# print "A normal", A
		# print "A inv -", np.linalg.inv(A)
		x[i] = np.dot(np.linalg.inv(A) , b)

	return x #returns the real points


def condition(P1,P2,DistA,DistB):
    # print (numpy.linalg.norm(P2-P1)- DistA)
    # print " < "
    # print DistB
    # print " < "
    # print (numpy.linalg.norm(P2-P1) + DistA)
    # print "////////////"
    if( ((np.linalg.norm(P2-P1)- DistA) < DistB)
        and (DistB < (np.linalg.norm(P2-P1) + DistA))):
        	print "OK"
		return True
    else:
        	print "NOT OK"
		return False

def full_condition(P1,P2,P3,DistA,DistB,DistC):
	A = condition(P1,P2,DistA,DistB)
	B = condition(P1,P3,DistA,DistC)
	C = condition(P2,P3,DistB,DistC)

	return A and B and C

def km_dist(p, lat, lon):
	return mpu.haversine_distance((LOC_mtz[p][0], LOC_mtz[p][1]) , (lat,lon))

def tri_lat(p,i,j,k):
	convert_np() #initializate BTS_mtz
	dist_mtz = dist_array() #m
	

	P1 = np.copy(BTS_mtz[i])
	P2 = np.copy(BTS_mtz[j])
	P3 = np.copy(BTS_mtz[k])
	print 

	#Translade circles references (BTS)
	ex = (P2 - P1)/(np.linalg.norm(P2 - P1))
	i = np.dot(ex,P3 - P1)
	ey = (P3 - P1 - i*ex)/(np.linalg.norm(P3 - P1 - i*ex))
	ez = np.cross(ex,ey)
	d = np.linalg.norm(P2 - P1)
	j = np.dot(ey,P3 - P1)

	Dist1 = np.copy(dist_mtz[int(p)][int(i)])
	Dist2 = np.copy(dist_mtz[int(p)][int(j)])
	Dist3 = np.copy(dist_mtz[int(p)][int(k)])

	print np.array([Dist1,Dist2,Dist3])

	print full_condition(P1,P2,P3,Dist1,Dist2,Dist3)
	#plug and chug
	
	x = (pow(Dist1,2) - pow(Dist2,2) + pow(d,2))/(2*d)
	y = (( pow(Dist1,2) - pow(Dist3,2) + pow(i,2) + pow(j,2))/(2*j))
	- ((i/j) *x)

	#only contendoa
	test = (pow(Dist1,2) - pow(x,2) - pow(y,2))
	if(test < 0):
		z = 1
	else:
		z = np.sqrt(test)

	

	#triPT dist_array
	triPt = P1 + x*ex + y*ey + z*ez
	
	#Back
	earthR = 6371

	lat = math.degrees(math.asin(triPt[2] / earthR))
	lon = math.degrees(math.atan2(triPt[1],triPt[0]))

	print LOC_mtz[p][0], LOC_mtz[p][1]
	print lat, lon

	dist = km_dist(p,lat,lon)
	print dist
		


if __name__ == "__main__":
	# try1 = rss_lat(0,1,2,3)
	#print "pred0 -",try1[0]
	print "test"
	#TODO: check if the point has 3 circles overlapping. if so, the best comb?
	#TODO: write .csv file base https://github.com/BrianSanderson/trilateration/blob/master/trilat.py

	#final point, bstindexs
	tri_lat(2,0,1,2)
