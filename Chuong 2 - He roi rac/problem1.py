#.....................................................................
# code Python cho phương pháp phần tử hữu hạn 
# problem1.py

# thêm các thư viện cần thiết
import sys	
sys.path.insert(0,'E:/My file/CODE FEM/Source code/Function')

import numpy as np
from outputDisplacementReactions import*

# các nút phần tử: được kết nối tại các phần tử
elementNodes=np.matrix('1 2;2 3;2 4')

# numberElements: số phần tử
numberElements=np.size(elementNodes, axis=(0))

# numberNodes: số nút
numberNodes=int(4)

# cho cấu trúc: 
	# displacements: vecto chuyển vị 
	# force : vecto lực
	# độ cứng: ma trận độ cứng
displacements=np.zeros((numberNodes,1))
force=np.zeros((numberNodes,1))
stiffness=np.zeros((numberNodes,numberNodes))

# Cho lực tại nút 2
nodeforce=2
force[nodeforce-1]=10.0

# Tính toán ma trận độ cứng toàn hệ

for i in range (0,numberElements) :
	# elementDof: bậc tự do (Dof)
	elementDof=elementNodes[i]
	idx=np.ix_([elementDof[0,0]-1,elementDof[0,1]-1],[elementDof[0,0]-1,elementDof[0,1]-1])
	stiffness[idx]=stiffness[idx]+[[1,-1],[-1,1]]
	
# điều kiện biên và cách giải 
# Tổng số bậc tự do 
GDof=4
GDofMatrix=[]
for i in range (0,GDof) :
	GDofMatrix.append([i])

# quy định bậc tự do bị ràng buộc
prescribedDof=np.matrix('1;3;4')

# free Dof : activeDof 
activeDof=np.setdiff1d(GDofMatrix,(prescribedDof-1))

# cách giải 
displacements=np.divide(force[activeDof],stiffness[activeDof,activeDof])

# định vị tất cả các chuyển vị
displacements1=np.zeros((numberNodes,1))
displacements1[activeDof]=displacements

# chuyển vị/lực đầu ra
outputDisplacementReactions(displacements1,stiffness,GDof,prescribedDof)