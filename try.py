#!/usr/bin/python
#Similar to TOA
#Try to write a program to solve the exercise of Trilateration GO

import csv
import numpy as np
import math

#Constants
data_BST = csv.reader(open('Dados/dados_BTSs.csv', 'rb'), delimiter=",", quotechar='|')
data_LOC = csv.reader(open('Dados/LocTreino_Equipe_4.csv', 'rb'))
#1501 linhas -> 1500 pontos

array_lon, array_lat = [], []
lat_index = 1
lon_index = 2
d0 = 1	#must be greater than 0, ideal 1, distance target-target
B = 6 #number of BSTs
L = 2 #Dimension 2D

#TODO: Decisions and completness test
def read_csv():
	#nem precisa dos arrays, da pra ir direto no BST
	for row in data_BST:
		array_lat.append(row[lat_index])	
		array_lon.append(row[lon_index])
	array_lat.pop(0)
	array_lon.pop(0)

	#print array_lat
	#print array_lon

def convert_np():
	#pode ser uma tupla
	BST_mtz = np.zeros([B,L])
	for i in range(0,B):
		BST_mtz[i][0] = array_lat[i]	
		BST_mtz[i][1] = array_lon[i]
	#print BST_mtz
	return BST_mtz

#Array of distances for each point
def estimate_dst():
	loc_list = list(data_LOC)
	loc_list.pop(0);

	LOC_mtz = np.zeros([len(loc_list),L])	
	dist_mtz = np.zeros([len(loc_list),B])			
	
	for i in range(0,len(loc_list)):
		LOC_mtz[i][0] = loc_list[i][1]
		LOC_mtz[i][1] = loc_list[i][2] 
		for j in range (0,6):
			temp = loc_list[i][j+3]
			dist_mtz[i][j]  = 1*pow(10, (float(temp)/(10*B)))
			#dist_mtz[i][j] = temp

	
	print LOC_mtz[0]
	print dist_mtz[0]
	return dist_mtz
	
def points(i_index,j_index,k_index, p_index):
	BST_mtz = convert_np()
	#print BST_mtz
	mtz_aux = np.zeros([L,L])

	mtz_aux[0][0] = BST_mtz[j_index][0] - BST_mtz[i_index][0]
	mtz_aux[0][1] = BST_mtz[j_index][1] - BST_mtz[i_index][1]
	mtz_aux[1][0] = BST_mtz[k_index][0] - BST_mtz[i_index][0]
	mtz_aux[1][1] = BST_mtz[k_index][1] - BST_mtz[i_index][1]
	
	dist_mtz = estimate_dst()
	array_aux = np.array([0,0])
	temp1 = pow(np.linalg.norm(BST_mtz[j_index]),2)
	temp2 = pow(np.linalg.norm(BST_mtz[i_index]),2)
	temp3 = pow(dist_mtz[p_index][j_index],2) - pow(dist_mtz[p_index][i_index],2)

	array_aux[0] = temp1 - temp2 - temp3
	
	temp1 = pow(np.linalg.norm(BST_mtz[k_index]),2)
	temp3 = pow(dist_mtz[p_index][k_index],2) - pow(dist_mtz[p_index][i_index],2)

	
	array_aux[1] = temp1 - temp2 - temp3
	
	
 	#array_aux[1] = 
	print array_aux
	
	#print mtz_aux


						
if __name__ == "__main__":
	read_csv()
	points(0,1,2,1)	
