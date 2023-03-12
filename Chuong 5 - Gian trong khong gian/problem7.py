#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem7.py

# thêm các thư viện cần thiết

import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code')

import numpy as np
import math as m
import matplotlib.pyplot as plt
from solution import *
from outputDisplacementReactions import *
from formStiffness3Dtruss import *
from stresses3Dtruss import *

# E : mô đun đàn hồi
# A: diện tích mặt cắt ngang
E=1.2E6
A=np.matrix('0.302;0.729;0.187')

# các hệ tọa độ và liên kết
elementNodes=np.matrix('1 2;1 3;1 4')
nodeCoordinates=np.matrix('72 0 0; 0 36 0; 0 36 72; 0 0 -48')
numberElements=np.size(elementNodes,0)
numberNodes=np.size(nodeCoordinates,0)
numberspring=0
xx=np.zeros((1,numberNodes-numberspring))
yy=np.zeros((1,numberNodes-numberspring))
zz=np.zeros((1,numberNodes-numberspring))
for i in range(0,numberNodes-numberspring):
	xx[0,i]=nodeCoordinates[i,0]
	yy[0,i]=nodeCoordinates[i,1]
	zz[0,i]=nodeCoordinates[i,2]

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
GDof=3*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
force=np.zeros((GDof,1))

# đặt tải trọng tại bậc tự do thứ 3
Dofforce=3
force[Dofforce-1]=-1000

# Ma trận độ cứng 
stiffness=formStiffness3Dtruss(GDof,numberElements,elementNodes,xx,yy,zz,E,A)
	
# điều kiện biên và cách giải
prescribedDof=[[2]]
for i in range(4,13):
	prescribedDof.append([i])
prescribedDof=np.matrix(prescribedDof)

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# vẽ các chuyển vị

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

# Ứng suất ở mỗi phần tử 
stresses3Dtruss(numberElements,elementNodes,xx,yy,zz,displacements,E)