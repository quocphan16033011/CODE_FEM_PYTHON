import numpy as np
import math 
def formStiffnessBernoulliBeam(GDof,numberElements,elementNodes,numberNodes,xx,EI,P):
	force=np.zeros((GDof,1))
	stiffness=np.zeros((GDof,GDof))
	for i in range(0,numberElements) :
		indice=elementNodes[i]
		elementDof=[2*(indice[0]-1),2*(indice[1]-1)-1,2*(indice[1]-1),2*(indice[1]-1)+1]
		length_element=xx[indice[1]-1,0]-xx[indice[0]-1,0]
		k1=EI/(length_element**3)*np.matrix([[12,6*length_element,-12,6*length_element],[6*length_element,4*length_element**2,-6*length_element,2*length_element**2],[-12,-6*length_element,12,-6*length_element],[6*length_element,2*length_element**2,-6*length_element,4*length_element**2]])
		f1=np.matrix([P*length_element/2,P*length_element**2/12,P*length_element/2,-P*length_element**2/12]).T
		idx_1=np.ix_(elementDof)
		force[idx_1]=force[idx_1]+f1
		idx=np.ix_(elementDof,elementDof)
		stiffness[idx]=stiffness[idx]+k1
	return stiffness,force
