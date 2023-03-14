#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem4.py

# thêm các thư viện cần thiết
import sys
sys.path.insert(0,'E:/My file/CODE FEM/Source code/Function')

import numpy as np
import matplotlib.pyplot as plt
from solution import *
from outputDisplacementReactions import *
from formStiffness2Dtruss import *
from stresses2Dtruss import *

# E : mô đun đàn hồi
# A: diện tích mặt cắt ngang
E=30E6
A=2
EA=E*A

# các hệ tọa độ và liên kết
numberElements=3
numberNodes=4
elementNodes=np.matrix('1,2;1,3;1,4')
nodeCoordinates=np.matrix('0,0;0,120;120,120;120,0')
xx=np.zeros((1,numberNodes))
yy=np.zeros((1,numberNodes))
for i in range(0,numberNodes):
	xx[0,i]=nodeCoordinates[i,0]
	yy[0,i]=nodeCoordinates[i,1]

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
GDof=2*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
force=np.zeros((GDof,1))

# đặt tải trọng bậc tự do thứ 2
Dofforce=2
force[Dofforce-1]=-10000

# Ma trận độ cứng 
stiffness=formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA)

# điều kiện biên và cách giải
prescribedDof=[]
for i in range(3,9):
	prescribedDof.append([i])
prescribedDof=np.matrix(prescribedDof)

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# vẽ các chuyển vị

# Ứng suất ở mỗi phần tử
stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E)

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

