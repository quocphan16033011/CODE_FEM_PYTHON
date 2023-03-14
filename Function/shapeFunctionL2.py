# Các hàm dạng và các đạo hàm liên quan đến tọa độ tự nhiên 
import numpy as np
def shapeFunctionL2(xi) :
	# hàm hình dạng và dẫn xuất cho các phần tử L2
	# shape : cách hàm hình dạng
	# naturalDerivatives: các dẫn xuất w.r.t. xi
	# xi: tọa độ tự nhiên (-1 ... +1)
	shape=[1-xi,1+xi]
	shape=(np.matrix(shape)/2).T
	naturalDerivaties=[[-1],[1]]
	naturalDerivaties=np.matrix(naturalDerivaties)/2
	return shape,naturalDerivaties
