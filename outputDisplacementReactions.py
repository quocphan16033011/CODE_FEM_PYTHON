# Hàm xuất ra kết quả chuyển vị và phản lực liên kết của thanh thẳng chịu kéo nén đúng tâm

def outputDisplacementReactions(displacements,stiffness,GDof,prescribedDof) :
	# xóa bộ nhớ
	import sys
	sys.modules[__name__].__dict__.clear()
	#Thêm các thư viện 
	import numpy as np	
	force=np.matmul(stiffness,displacements)
	reaction=force[prescribedDof-1]
	print(end="\n")
	print('Chuyển vị')
	print(end="\n")
	print('','Bậc tự do','    ','Giá trị')
	for i in range(0,len(displacements)):
		if i <9 :
			print('   ',i+1,'          ',displacements[i])
		elif i<99 : print('   ',i+1,'         ',displacements[i])
		else : print('   ',i+1,'        ',displacements[i])
	print(end="\n")
	print('Phản lực liên kết',end="\n")
	print(end="\n")
	print('','Bậc tự do','    ','Giá trị')
	for i in range(0,len(prescribedDof)):
		if prescribedDof[i,0] <10 :
			print('   ',prescribedDof[i,0],'          ',reaction[i,0])
		elif prescribedDof[1,0] < 99 : print('   ',prescribedDof[i,0],'         ',reaction[i,0])
		else :  print('   ',prescribedDof[i,0],'        ',reaction[i,0])
	print(end="\n")