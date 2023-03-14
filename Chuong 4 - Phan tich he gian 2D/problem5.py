#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem5.py

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
E=70000
A=300
EA=E*A

# các hệ tọa độ và liên kết
numberElements=11
numberNodes=6
elementNodes=np.matrix('1 2;1 3;2 3;2 4;1 4;3 4;3 6;4 5;4 6;3 5;5 6')
nodeCoordinates=np.matrix('0 0;0 3000;3000 0;3000 3000;6000 0;6000 3000')
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

# đặt tải trọng tại bậc tự do thứ 4,8,12
Dofforce_1=4
Dofforce_2=8
Dofforce_3=12
force[Dofforce_1-1]=-50000
force[Dofforce_2-1]=-100000
force[Dofforce_3-1]=-50000

# Ma trận độ cứng 
stiffness=formStiffness2Dtruss(GDof,numberElements,elementNodes,numberNodes,nodeCoordinates,xx,yy,EA)

# điều kiện biên và cách giải
prescribedDof=np.matrix('1;2;10')

# cách giải
displacements=solution(GDof,prescribedDof,stiffness,force)

# vẽ các chuyển vị

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof)

# Ứng suất ở mỗi phần tử 
stresses2Dtruss(numberElements,elementNodes,xx,yy,displacements,E)