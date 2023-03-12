#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem6.py

# thêm các thư viện cần thiết
import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code')

import numpy as np
import math as m
import matplotlib.pyplot as plt
from solution import *
from outputDisplacementReactions import *
from formStiffness2Dtruss import *
from stresses2Dtruss import *

# E : mô đun đàn hồi
# A: diện tích mặt cắt ngang
E=210000
A=500
EA=E*A

# các hệ tọa độ và liên kết
elementNodes=np.matrix('1 2;1 3;1 4')
nodeCoordinates=np.matrix([[0,0],[-5000*m.cos(m.pi/4),5000*m.sin(m.pi/4)],[-10000,0]])
numberElements=np.size(elementNodes,0)
numberNodes=np.size(elementNodes,0)+1
numberspring=1
xx=np.zeros((1,numberNodes-numberspring))
yy=np.zeros((1,numberNodes-numberspring))
for i in range(0,numberNodes-numberspring):
	xx[0,i]=nodeCoordinates[i,0]
	yy[0,i]=nodeCoordinates[i,1]

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng

GDof=2*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
force=np.zeros((GDof,1))
stiffness=np.zeros((GDof,GDof))

# đặt tải trọng tại bậc tự do thứ 2
Dofforce=2
force[Dofforce-1]=-25000

# Ma trận độ cứng 
for i in range(0,numberElements-numberspring) :
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
idx1=np.ix_([1,6],[1,6])
stiffness[idx1]=stiffness[idx1]+ np.dot(2000,[[1,-1],[-1,1]])

# điều kiện biên và cách giải
prescribedDof=[]
for i in range(3,9):
	prescribedDof.append([i])
prescribedDof=np.matrix(prescribedDof)

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# vẽ các chuyển vị

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

# Ứng suất ở mỗi phần tử 
sigma=np.zeros((1,numberElements))
for i in range(0,numberElements-numberspring) :
	indice=elementNodes[i]
	elementDof=[indice[0,0]*2-2, indice[0,0]*2-1, indice[0,1]*2-2, indice[0,1]*2-1]
	xa=xx[0,indice[0,1]-1]-xx[0,indice[0,0]-1]
	ya=yy[0,indice[0,1]-1]-yy[0,indice[0,0]-1]
	length_element=math.sqrt(xa*xa+ya*ya)
	c=xa/length_element
	s=ya/length_element
	idx=np.ix_(elementDof)
	sigma[0,i]=E/length_element*np.matrix([-c,-s,c,s])*displacements[idx]	
print('Ứng suất')
print(end='\n')
print('','Phần tử','      ','Giá trị')
for i in range(0,numberElements-numberspring):
	print('   ',i+1,'          ',sigma.T[i])
print(end='\n')