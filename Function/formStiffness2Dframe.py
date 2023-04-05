import numpy as np
import math as m

def formStiffness2Dframe(GDof,numberElements,elementNodes,numberNodes,xx,yy,EI,EA):
	stiffness=np.zeros((GDof,GDof))
	for i in range(0,numberElements):
		indice=elementNodes[i]
		elementDof=np.hstack((indice-1,indice+numberNodes-1,indice+2*numberNodes-1))
		nn=len(indice)
		xa=xx[indice[1]-1]-xx[indice[0]-1]
		ya=yy[indice[1]-1]-yy[indice[0]-1]
		length_element=m.sqrt(xa**2+ya**2)
		cosa=xa/length_element
		sina=ya/length_element
		L=np.vstack((np.hstack((cosa*np.eye(2),sina*np.eye(2),np.zeros((2,2)))),np.hstack((-sina*np.eye(2),cosa*np.eye(2),np.zeros((2,2)))),np.hstack((np.zeros((2,4)),np.eye(2)))))
		oneu=np.matrix([[1,-1],[-1,1]])
		oneu2=np.matrix([[1,-1],[1,-1]])
		oneu3=np.matrix([[1,1],[-1,-1]])
		oneu4=np.matrix([[4,2],[2,4]])
		k1=np.vstack((np.hstack((EA/length_element*oneu,np.zeros((2,4)))),np.hstack((np.zeros((2,2)),12*EI/length_element**3*oneu,6*EI/length_element**2*oneu3)),np.hstack((np.zeros((2,2)),6*EI/length_element**2*oneu2,EI/length_element*oneu4))))
		idx=np.ix_(elementDof,elementDof)
		stiffness[idx]=stiffness[idx]+np.dot(np.dot(L.T,k1),L)

	return stiffness