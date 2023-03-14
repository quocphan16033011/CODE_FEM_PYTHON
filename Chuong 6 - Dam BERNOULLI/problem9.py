#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem7.py

# thêm các thư viện cần thiết

import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code/Function')

import numpy as np
import math as m
import matplotlib.pyplot as plt
from solution import *
from outputDisplacementReactions import *
from formStiffnessBernoulliBeam import *

# E : mô đun đàn hồi
# I: Moment thứ 2 của khu vực
E=1
I=1
EI=E*I

# các hệ tọa độ và liên kết
numberElements=80
nodeCoordinates=np.transpose(np.matrix(np.linspace(0,1,numberElements+1)))
numberNodes=np.size(nodeCoordinates,0)
numberspring=0
xx=np.zeros((numberNodes-numberspring,1))
for i in range(0,numberNodes-numberspring):
	xx[i]=nodeCoordinates[i]
elementNodes=np.zeros((numberElements,2),dtype=int)
for i in range(0,numberElements) :
	elementNodes[i,0]=i+1
	elementNodes[i,1]=i+2

# Lực phân bố
P=-1

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
GDof=2*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
force=np.zeros((GDof,1))
stiffness=np.zeros((GDof,GDof))

# Ma trận độ cứng và tải trọng tác dụng lên dầm 
[stiffness,force]=formStiffnessBernoulliBeam(GDof,numberElements,elementNodes,numberNodes,xx,EI,P)

# điều kiện biên và cách giải
fixedNodeU=np.matrix([1,2*numberElements+1],dtype=int).T
prescribedDof=fixedNodeU

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# Vẽ hình dạng hình học
U=[]
for i in range (0,GDof) :
	if i % 2 ==0 :		
		U.append(displacements[i])
plt.scatter([nodeCoordinates[:,0]],U,s=3.5)
plt.show()

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)