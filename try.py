#!/usr/bin/python
#Similar to TOA
#Try to write a program to solve the exercise of Trilateration GO

#useful http://www.meridianoutpost.com/resources/etools/calculators/calculator-latitude-longitude-distance.php?
#http://www.ambrsoft.com/TrigoCalc/Circles3/Intersection.htm
#https://gis.stackexchange.com/questions/66/trilateration-using-3-latitude-longitude-points-and-3-distances/415#415
#https://en.wikipedia.org/wiki/Trilateration
#http://obeattie.github.io/gmaps-radius/?lat=51.500358&lng=-0.125506&z=10&u=mi&r=5

import csv
import numpy as np
import math
import mpu
#pip install mpu --user
#pip install enum --user
# import pyproj as pj

from PyRadioLoc.Pathloss.Models import FreeSpaceModel
from PyRadioLoc.Pathloss.Models import FlatEarthModel
# from PyRadioLoc.Pathloss.Models import LeeModel
# from PyRadioLoc.Pathloss.Models import EricssonModel
from PyRadioLoc.Pathloss.Models import Cost231Model
from PyRadioLoc.Pathloss.Models import Cost231HataModel
# from PyRadioLoc.Pathloss.Models import OkumuraHataModel
# from PyRadioLoc.Pathloss.Models import Ecc33Model
from PyRadioLoc.Pathloss.Models import SuiModel


#Data
data_BTS = csv.reader(open('Dados/dados_BTSs.csv', 'rb'))
# data_LOC = csv.reader(open('Dados/LocTreino_Equipe_4.csv', 'rb'))
data_LOC = csv.reader(open('Dados/LocTest.csv', 'rb'))

#Constants
lat_index = 1
lon_index = 2
d0 = 1	#must be greater than 0
PL = 2 #Pl expoent in free space
factor = 8 #For FPSL aproximate
B = 6 #number of BTSs
L = 3 #Dimension 2D
P = 1500 #number of points
TestP = 200 #points for test
GSM = 1800 #frequency for FPSL MHZ and meters (-27.55)https://en.wikipedia.org/wiki/Free-space_path_loss

# LOC_mtz = np.zeros([P,L]) #Real points
LOC_mtz = np.array([[1.]*L]*TestP)
# BTS_mtz = np.zeros([B,L]) #Bts Points
BTS_mtz = np.array([[1.]*L]*B)

#results
#pontoId, lat, lon , lat_pred, lon_pred e erro_loc. (metros)
RES_list = []
RES_list.append(['pontoId', 'lat','lon','lat_pred','lon_pred','erro_loc'])

list_bts = list(data_BTS)
list_bts.pop(0)

#TODO: Decisions and completness test
def coord_conv(lat,lon,alt):
	lat = float(lat)
	lon = float(lon)
	earthR = 6371
	x = earthR *(math.cos(math.radians(lat)) * math.cos(math.radians(lon)))
	y = earthR *(math.cos(math.radians(lat)) * math.sin(math.radians(lon)))
	z = earthR *(math.sin(math.radians(lat)))

	return np.array([x,y,z])

def km_dist(p, lat, lon):
	return mpu.haversine_distance((LOC_mtz[p][0], LOC_mtz[p][1]) , (lat,lon))

def convert_np():
	# list_bts = list(data_BTS)
	# list_bts.pop(0)

	for i in range(0,B):
		# BTS_mtz[i] = gps_ecef(list_bts[i][1], list_bts[i][2], 0)
		print "number: ", i, "-", list_bts[i][1],",",list_bts[i][2]
		BTS_mtz[i] = coord_conv(list_bts[i][1], list_bts[i][2], 0)
		# for j in range (0,L):
		# 	BTS_mtz[i][j] = list_bts[i][j+1]

	# print BTS_mtz


def ML_estimate(vi):
	return d0*pow(10 , (float(vi)/(10*PL)))

def FPSL_estimate(vi):
	#FSPL(db), d in meters, f in MHZ , cte = -27.55
	#FSPL(db), d in meters, f in GHZ , cte = -87.55
	#FSPL(db), d in Kmeters, f in MHZ , cte = +32.45
	#FSPL(db), d in Kmeters, f in GHZ , cte = +92.45
	f = 1800 #original MHZ
	cte = +32.55
	#double sqrt seems to be the best one
	return pow((pow(10, (( float(vi)-cte)-(20*math.log(f,10)))/20)), 1./(factor*2))

def SUI_estimate(vi):
	m1 = SuiModel(1901)
	return m1.dist_path(float(vi))

def Cost231Hata_estimate(vi):
	m1 = Cost231HataModel(1800)
	return m1.dist_path(float(vi))

#return distance to target array and modify the real points
diff = np.array([[1.]*B]*TestP)
def dist_array():

	#Read .csv file with 1500 points
	loc_list = list(data_LOC)
	loc_list.pop(0)
	print loc_list[0]

	dist_mtz = np.array([[1.]*B]*TestP)

	#Each point will a array of 6 distance (anchor - target)
	for i in range(0,TestP):
		# LOC_mtz[i] = gps_ecef(loc_list[i][1], loc_list[i][2],0) #lat, lon ,hei
		RES_list.append([(loc_list[i][0])])
		#LOC_mtz[i] = coord_conv(loc_list[i][1], loc_list[i][2],0)
		for j in range (0,2):
			LOC_mtz[i][j] = float(loc_list[i][j+1])
		for j in range (0,B):
			diff[i][j] = km_dist(i, float(list_bts[j][1]), float(list_bts[j][2]))
			# print loc_list[i][j+3]
			if(loc_list[i][j+3] != 'NA'):
				#dist_mtz[i][j] = Cost231Hata_estimate(loc_list[i][j+3])
				#dist_mtz[i][j] = SUI_estimate(loc_list[i][j+3])
				dist_mtz[i][j] = FPSL_estimate(loc_list[i][j+3])
				#dist_mtz[i][j] = (Cost231Hata_estimate(loc_list[i][j+3]) + SUI_estimate(loc_list[i][j+3]))/2

	#print "Dist ponto 1 -",  dist_mtz[0]
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

#TODO: intersect in 2 points. best case, each
def condition(P1,P2,DistA,DistB):
    # print (numpy.linalg.norm(P2-P1)- DistA)
    # print " < "
    # print DistB
    # print " < "
    # print (numpy.linalg.norm(P2-P1) + DistA)
    # print "////////////"
	if((DistA - DistB) < np.linalg.norm(P2-P1) and (DistA + DistB) > np.linalg.norm(P2-P1)):
		return True
	else:
		return False
    # if( ((np.linalg.norm(P2-P1)- DistA) < DistB)
    #     and (DistB < (np.linalg.norm(P2-P1) + DistA))):
    #     	print "OK"
	# 	return True
    # else:
    #     	# print "NOT OK"
	# 	return False

def full_condition(P1,P2,P3,DistA,DistB,DistC):
	A = condition(P1,P2,DistA,DistB)
	B = condition(P1,P3,DistA,DistC)
	C = condition(P2,P3,DistB,DistC)

	return A and B and C

def min_radius(P1,P2,P3,D1,D2,D3):
	if(min(D1,D2,D3) == D1):
		return P1[2]
	elif(min(D2,D2,D3) == D2):
		return P2[2]
	else:
		return P3[2]

def tri_lat(p_index,i_index,j_index,k_index, dist_mtz):
	#print p_index, i_index, j_index, k_index
	# convert_np() #initializate BTS_mtz
	# dist_mtz = dist_array() #m


	P1 = np.copy(BTS_mtz[i_index])
	P2 = np.copy(BTS_mtz[j_index])
	P3 = np.copy(BTS_mtz[k_index])

	Dist1 = np.copy(dist_mtz[p_index][i_index])
	Dist2 = np.copy(dist_mtz[p_index][j_index])
	Dist3 = np.copy(dist_mtz[p_index][k_index])

	#Translade circles references (BTS)
	ex = (P2 - P1)/(np.linalg.norm(P2 - P1))
	i = np.dot(ex,P3 - P1)
	ey = (P3 - P1 - i*ex)/(np.linalg.norm(P3 - P1 - i*ex))
	ez = np.cross(ex,ey)
	d = np.linalg.norm(P2 - P1)
	j = np.dot(ey,P3 - P1)


	#plug and chug
	x = (pow(Dist1,2) - pow(Dist2,2) + pow(d,2))/(2*d)
	y = (( pow(Dist1,2) - pow(Dist3,2) + pow(i,2) + pow(j,2))/(2*j))
	- ((i/j) *x)

	#Why y got so penalized i,j,k

	#only contendoa
	test = (pow(Dist1,2) - pow(x,2) - pow(y,2))
	if(test < 0):
		z = min(x,y)
		#z = 1
		# print "Forced Z"
		# return float("inf")
		return (float("inf"),0,0)
	else:
		z = np.sqrt(test)

	# print "Cartesian Pred - ", np.array([x,y,z]), "\n"
	#triPT dist_array
	triPt = P1 + x*ex + y*ey + z*ez

	#Back
	earthR = 6371
	if(triPt[2]/earthR < -1 or triPt[2]/earthR > 1):
		# print "Forced again"
		# return float("inf")
		return (float("inf"),0,0)
	lat = math.degrees(math.asin(triPt[2] / earthR))
	lon = math.degrees(math.atan2(triPt[1],triPt[0]))

	# print "REal - ",LOC_mtz[p_index][0], LOC_mtz[p_index][1]
	# print "Pred - ",lat, lon

	dist = km_dist(p_index,lat,lon)
	# print "Erro -", dist*1000, "m"
	return (dist,lat,lon)


#Menor erro com determinada BTS
def min_error(p_index, bst_index, dist_mtz):
	list = []
	for i in range (0,6):
		if(i!=bst_index):
			for j in range (0,6):
				if(j!=bst_index and j!= i):
					list.append( tri_lat(p_index,bst_index,i,j,dist_mtz))
					# print "i,j,k -", bst_index,i,j
					#print "C,I,J: ",bst_index, i,j
	#print "Min Dist - ", min(list)*1000, "m"
	return min(list)

#Menor erro com combinacao de todas as BTSs
def best_error(p_index,dist_mtz):
	list = []
	for i in range (0,6):
		list.append(min_error(p_index,i,dist_mtz))

	return min(list)

def all_points(dist_mtz, points):
	sum = 0
	p = points
	for i in range(0,TestP):
		(err,lat,lon) = best_error(i,dist_mtz)
		RES_list[i+1].append(LOC_mtz[i][0])
		RES_list[i+1].append(LOC_mtz[i][1])
		RES_list[i+1].append(lat)
		RES_list[i+1].append(lon)
		RES_list[i+1].append(err*1000)
		# print "Point:", i, " Min Error ", err*1000, "m"
		# print "Real: ", LOC_mtz[i][0],",",LOC_mtz[i][1]
		# print "Pred: ", lat,",",lon, "\n"
		if(err != float("inf")):
			sum += err
		else:
			p = p-1

	print "Average error - ", (sum/p)*1000, "m" , "with:", p ," points"
if __name__ == "__main__":
	# try1 = rss_lat(0,1,2,3)
	#print "pred0 -",try1[0]
	print "test"
	#TODO: check if the point has 3 circles overlapping. if so, the best comb?
	#TODO: write .csv file base https://github.com/BrianSanderson/trilateration/blob/master/trilat.py

	#final point, bstindexs INIT
	#TODO: Apenas em uma leitura.
	convert_np()
	dist_mtz = dist_array() #m

	#tri_lat(0,1,2,4, dist_mtz)
	print best_error(0,dist_mtz)[0]*1000,"m - Min error"
	print "R-",diff[0]
	print "P-",dist_mtz[0]
	print "Real: ", LOC_mtz[0][0],",",LOC_mtz[0][1]

	all_points(dist_mtz,200)
	#print RES_list

	# with open("Dados/Resultados_3.csv", "wb") as f:
	# 	writer = csv.writer(f)
	# 	writer.writerows(RES_list)
