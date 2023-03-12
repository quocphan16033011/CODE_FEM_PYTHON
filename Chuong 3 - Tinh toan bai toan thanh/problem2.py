#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem3.py

# thêm các thư viện cần thiết
import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code')

import numpy as np
from shapeFunctionL2 import *
from solution import *
from outputDisplacementReactions import *

# E : mô đun đàn hồi
# A: diện tích mặt cắt ngang
# L: chiều dài thanh
E=int(30E6)
A=int(1)
EA=E*A
L=int(90)

# các hệ tọa độ và liên kết

# numberElements: số phần tử
numberElements=int(3)
nodeCoordinates=np.linspace(0,L,numberElements+1)
xx=nodeCoordinates

# numberNodes: số nút
numberNodes=len(nodeCoordinates)

# elementNodes: các liên kết trong phần tử
ii=range(0,numberElements)
elementNodes=[]
for i in ii :
	elementNodes.append([i+1,i+2])

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
displacements=np.zeros((numberNodes,1))
force=np.zeros((numberNodes,1))
stiffness=np.zeros((numberNodes,numberNodes))

# Đặt tải trọng tại nút 2
nodeforce=2
force[nodeforce-1]=int(3000)

# Tính toán ma trận độ cứng toàn hệ
for i in range (0,numberElements) :
	# elementDof: bậc tự do (Dof)
	elementDof=elementNodes[i]
	nn=len(elementDof)
	length_element=nodeCoordinates[elementDof[1]-1]-nodeCoordinates[elementDof[0]-1]
	detJacobian=length_element/2
	invJacobian=1/detJacobian
	# điểm Gauss trung tâm (xi=0, weight W=2)
	[shape,naturalDerivaties]=shapeFunctionL2(0)
	Xderivatives=np.dot(naturalDerivaties,invJacobian)
	#ma trận B
	B=np.zeros((1,nn))
	for j in range(0,nn):
		B[0,j]=Xderivatives[j]
	idx=np.ix_([elementDof[0]-1,elementDof[1]-1],[elementDof[0]-1,elementDof[1]-1])
	stiffness[idx]=stiffness[idx]+np.dot(B.T,B)*2*detJacobian*EA

# điều kiện biên và cách giải
# bậc tự do quy định
fixedDof=[]
for i in range (0,len(nodeCoordinates)) :
	if nodeCoordinates[i]==np.min(nodeCoordinates) :
		fixedDof.append([i+1])
	elif nodeCoordinates[i]==np.max(nodeCoordinates) :
		fixedDof.append([i+1])

prescribedDof=np.matrix(fixedDof)

# Tổng số bậc tự do 
GDof=numberNodes

# chuyển vị của hệ
displacements=solution(GDof,prescribedDof,stiffness,force)

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

