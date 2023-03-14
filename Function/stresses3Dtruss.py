import numpy as np
import math
def stresses3Dtruss(numberElements,elementNodes,xx,yy,zz,displacements,E):
	sigma=np.zeros((1,numberElements))
	for i in range(0,numberElements) :
		indice=elementNodes[i]
		elementDof=[3*indice[0,0]-3, 3*indice[0,0]-2, 3*indice[0,0]-1, 3*indice[0,1]-3, 3*indice[0,1]-2, 3*indice[0,1]-1]
		x1=xx[0,indice[0,0]-1]
		y1=yy[0,indice[0,0]-1]
		z1=zz[0,indice[0,0]-1]
		x2=xx[0,indice[0,1]-1]
		y2=yy[0,indice[0,1]-1]
		z2=zz[0,indice[0,1]-1]
		L=math.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +(z2-z1)*(z2-z1))
		CXx = (x2-x1)/L
		CYx = (y2-y1)/L
		CZx = (z2-z1)/L
		idx=np.ix_(elementDof)
		sigma[0,i]=E/L*np.matrix([-CXx,-CYx,-CZx,CXx,CYx,CZx])*displacements[idx]	
	print('Ứng suất')
	print(end='\n')
	print('','Phần tử','      ','Giá trị')
	for i in range(0,numberElements):
		if i <9 :
			print('   ',i+1,'          ',sigma.T[i])
		else : print('   ',i+1,'         ',sigma.T[i])
	print(end='\n') 