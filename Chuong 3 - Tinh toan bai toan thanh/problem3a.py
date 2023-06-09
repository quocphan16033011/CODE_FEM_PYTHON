#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem3.py

# thêm các thư viện cần thiết
import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code/Function')

import numpy as np
from shapeFunctionL2 import *
from solution import *
from outputDisplacementReactions import *

# E : mô đun đàn hồi
# A : diện tích mặt cắt ngang
# k : độ cứng lò xo
E=int(7e4)
A=int(200)
EA=E*A
k=int(2000)

# các hệ tọa độ và liên kết

# numberElements: số phần tử
numberElements=int(3)
# numberNodes: số nút
numberNodes=int(4)
elementNodes=np.matrix('1 2;2 3;3 4')
nodeCoordinates=np.matrix('0 2000 4000 4000')
xx=nodeCoordinates

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
displacements=np.zeros((numberNodes,1))
force=np.zeros((numberNodes,1))
stiffness=np.zeros((numberNodes,numberNodes))

# Đặt tải trọng tại nút 2
nodeforce=2
force[nodeforce-1]=int(8000)

# Tính toán ma trận độ cứng toàn hệ
for i in range (0,numberElements) :
	# elementDof: bậc tự do (Dof)
	elementDof=elementNodes[i]
	nn=np.size(elementDof)
	if i<2 :
		length_element=nodeCoordinates[0,elementDof[0,1]-1]-nodeCoordinates[0,elementDof[0,0]-1]
		detJacobian=length_element/2
		invJacobian=1/detJacobian
		# điểm Gauss trung tâm (xi=0, weight W=2)
		[shape,naturalDerivaties]=shapeFunctionL2(0)
		Xderivatives=np.dot(naturalDerivaties,invJacobian)
		#ma trận B
		B=np.zeros((1,nn))
		for j in range(0,nn):
			B[0,j]=Xderivatives[j]
		idx=np.ix_([elementDof[0,0]-1,elementDof[0,1]-1],[elementDof[0,0]-1,elementDof[0,1]-1])
		stiffness[idx]=stiffness[idx]+np.dot(B.T,B)*2*detJacobian*EA
	else :
		idx=np.ix_([elementDof[0,0]-1,elementDof[0,1]-1],[elementDof[0,0]-1,elementDof[0,1]-1])
		stiffness[idx]=stiffness[idx]+np.dot(k,[[1,-1],[-1,1]])

# điều kiện biên và cách giải
# bậc tự do quy định
prescribedDof=np.matrix('1;4')

# Tổng số bậc tự do 
GDof=numberNodes

# chuyển vị của hệ
displacements=solution(GDof,prescribedDof,stiffness,force)

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

