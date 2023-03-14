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
# I : Moment thứ 2 của khu vực
# L : Chiều dài thanh
# k : Độ cứng lò xo
E=1E6
L=10
t=L/1000
I=1*t**3/12	
EI=E*I
k=10

# các hệ tọa độ và liên kết
numberElements=3
nodeCoordinates=np.transpose(np.matrix(np.linspace(0,L,numberElements+1)))
numberNodes=np.size(nodeCoordinates,0)
xx=np.zeros((numberNodes,1))
for i in range(0,numberNodes):
	xx[i]=nodeCoordinates[i]
elementNodes=np.zeros((numberElements,2),dtype=int)
for i in range(0,numberElements) :
	elementNodes[i,0]=i+1
	elementNodes[i,1]=i+2

# Lực phân bố
P=-1000

# cho cấu trúc:
	# displacements: vecto chuyển vị
	# force : vecto lực
	# stiffness: ma trận độ cứng
GDof=2*numberNodes #GDof: tổng số bậc tự do
displacements=np.zeros((GDof,1))
stiffnessSpring=np.zeros((GDof+1,GDof+1))
forceSpring=np.zeros((GDof+1,1))

# Ma trận độ cứng và tải trọng tác dụng lên dầm 
[stiffness,force]=formStiffnessBernoulliBeam(GDof,numberElements,elementNodes,numberNodes,xx,EI,P)

# Thêm lò xo vào
idx=np.ix_(np.arange(0,GDof,1),np.arange(0,GDof,1))
stiffnessSpring[idx]=stiffness
idx_1=np.ix_(np.arange(0,GDof,1))
forceSpring[idx_1]=force

idx_2=np.ix_([GDof-2,GDof],[GDof-2,GDof])
stiffnessSpring[idx_2]=stiffnessSpring[idx_2]+np.matrix([[k,-k],[-k,k]])

# điều kiện biên và cách giải
fixedNodeU=np.matrix([1],dtype=int).T
fixedNodeV=np.matrix([2],dtype=int).T
prescribedDof=np.vstack((fixedNodeU,fixedNodeV,GDof+1))

# cách giải
displacements=solution(GDof+1,prescribedDof,stiffnessSpring,forceSpring)

# Vẽ hình dạng hình học
U=[]
for i in range (0,GDof) :
	if i % 2 ==0 :		
		U.append(displacements[i])
plt.scatter([nodeCoordinates[:,0]],U,s=3.5)
#plt.show()

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements,stiffnessSpring,GDof+1,prescribedDof)

# Phương pháp giải chính xác của Bathe (Solutions Manual of Procedures ...)
load=np.matrix([[L*P/3],[L*P/3],[L*P/6]])
K=EI/(L**3)*np.matrix([[189,-108,27],[-108,135,-54],[27,-54,27+k*L**3/EI]])
X=np.dot(np.linalg.inv(K),load)
print('Phương pháp giải chính xác của Bathe (Solutions Manual of Procedures ...)')
print(end='\n')
print(X)