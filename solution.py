import numpy as np

# Hàm trả các vector chuyển vị

def solution(GDof,prescribedDof,stiffness,force) :
	GDofMatrix=[]
	for i in range (0,GDof) :
		GDofMatrix.append([i])

	activeDof=np.setdiff1d(GDofMatrix,prescribedDof-1)+1
	
	idx1=np.ix_(activeDof-1,activeDof-1)
	U=np.dot(np.linalg.inv(stiffness[idx1]),force[activeDof-1])
	displacements=np.zeros((GDof,1))
	displacements[activeDof-1]=U
	return displacements