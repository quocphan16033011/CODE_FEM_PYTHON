import numpy as np
import math
def axialinternalforce(numberElements,elementNodes,xx,yy,displacements,E,A):
	internal=np.zeros((1,numberElements))
	for i in range(0,numberElements) :
		indice=elementNodes[i]
		elementDof=[indice[0,0]*2-2, indice[0,0]*2-1, indice[0,1]*2-2, indice[0,1]*2-1]
		xa=xx[0,indice[0,1]-1]-xx[0,indice[0,0]-1]
		ya=yy[0,indice[0,1]-1]-yy[0,indice[0,0]-1]
		length_element=math.sqrt(xa*xa+ya*ya)
		c=xa/length_element
		s=ya/length_element
		idx=np.ix_(elementDof)
		internal[0,i]=E*A/length_element*np.matrix([-c,-s,c,s])*displacements[idx]	
	print('Nội lực dọc trục')
	print(end='\n')
	print('','Phần tử','      ','Giá trị')
	for i in range(0,numberElements):
		if i <9 :
			print('   ',i+1,'          ',internal.T[i])
		else : print('   ',i+1,'         ',internal.T[i])
	print(end='\n') 