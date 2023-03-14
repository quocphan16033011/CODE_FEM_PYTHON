import numpy as np
import math 
def formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA) :

	stiffness=np.zeros((GDof,GDof))
	for i in range(0,numberElements) :
		indice=elementNodes[i]
		elementDof=[indice[0,0]*2-2, indice[0,0]*2-1, indice[0,1]*2-2, indice[0,1]*2-1]
		xa=xx[0,indice[0,1]-1]-xx[0,indice[0,0]-1]
		ya=yy[0,indice[0,1]-1]-yy[0,indice[0,0]-1]
		length_element=math.sqrt(xa*xa+ya*ya)
		c=xa/length_element
		s=ya/length_element
		k1=round(EA/length_element,-1)*np.matrix([[c*c,c*s,-c*c,-c*s],[c*s,s*s,-c*s,-s*s],[-c*c,-c*s,c*c,c*s],[-c*s,-s*s,c*s,s*s]])
		idx=np.ix_(elementDof,elementDof)
		stiffness[idx]=stiffness[idx]+k1
	
	return stiffness