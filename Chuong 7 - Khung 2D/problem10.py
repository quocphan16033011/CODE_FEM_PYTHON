#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem10.py

# thêm các thư viện cần thiết

import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code/Function')

import numpy as np
import math as m
import matplotlib.pyplot as plt
from solution import *
from outputDisplacementReactions import *
from formStiffness2Dframe import *

# E : mô đun đàn hồi
# A : Diện tích mặt cắt ngang
# I: Moment thứ 2 của khu vực
E=210000
A=100
I=2E8
EA=E*A
EI=E*I

#Tạo tọa độ và các liên kết
numberElements=3
p1=3000*(1+m.cos(m.pi/4))
nodeCoordinates=np.matrix([[0,3000],[3000,3000],[p1,0],[p1+3000,0]],dtype=int)
elementNodes=np.zeros((numberElements,2),dtype=int)
for i in range(0,numberElements) :
	elementNodes[i,0]=i+1
	elementNodes[i,1]=i+2
numberNodes=np.size(nodeCoordinates,0)
xx=np.zeros((numberNodes,1),dtype=int)
yy=np.zeros((numberNodes,1),dtype=int)
for i in range(0,numberNodes):
	xx[i]=nodeCoordinates[i,0]
	yy[i]=nodeCoordinates[i,1]

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
GDof=3*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
force=np.zeros((GDof,1))

stiffness=formStiffness2Dframe(GDof,numberElements,elementNodes,numberNodes,xx,yy,EI,EA)

# Vector lực
Dofforce_1=6
Dofforce_2=7
Dofforce_3=10
Dofforce_4=11
force[Dofforce_1-1]=-10000
force[Dofforce_2-1]=-10000
force[Dofforce_3-1]=-5E6
force[Dofforce_4-1]=5E6

# điều kiện biên và cách giải
prescribedDof=np.matrix("1;4;5;8;9;12")

# Cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

#
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

idx=np.ix_([0,2,4])
idx_1=np.ix_([1,3,5])

#plt.plot(xx,yy,dashes=[4,5])
#plt.plot(xx-reaction[idx],yy-reaction[idx_1])
#plt.show()