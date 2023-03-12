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
E=210000
A=np.matrix('100;100;100;100')

# các hệ tọa độ và liên kết
elementNodes=np.matrix('1 2;1 3;1 4;1 5')
nodeCoordinates=np.matrix('4000 4000 3000;0 4000 0;0 4000 6000;4000    0 3000;8000 -1000 1000')
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

# đặt tải trọng tại bậc tự do thứ 2
Dofforce=2
force[Dofforce-1]=-10000

# Ma trận độ cứng 
stiffness=formStiffness3Dtruss(GDof,numberElements,elementNodes,xx,yy,zz,E,A)
	
# điều kiện biên và cách giải
prescribedDof=[]
for i in range(4,16):
	prescribedDof.append([i])
prescribedDof=np.matrix(prescribedDof)

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# vẽ các chuyển vị

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

# Ứng suất ở mỗi phần tử 
stresses3Dtruss(numberElements,elementNodes,xx,yy,zz,displacements,E)