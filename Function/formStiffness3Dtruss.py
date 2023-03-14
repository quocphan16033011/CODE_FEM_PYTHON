import numpy as np
import math 
def formStiffness3Dtruss(GDof,numberElements,elementNodes,xx,yy,zz,E,A) :

	stiffness=np.zeros((GDof,GDof))
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
		T = np.matrix([[CXx*CXx,CXx*CYx,CXx*CZx],[CYx*CXx,CYx*CYx,CYx*CZx],[CZx*CXx,CZx*CYx,CZx*CZx]])
		idx=np.ix_(elementDof,elementDof)
		T=np.vstack((np.hstack((T,-T)),np.hstack((-T,T))))
		stiffness[idx]=stiffness[idx]+E*A[i,0]/L*T
	
	return stiffness