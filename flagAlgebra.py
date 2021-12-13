import sys
import mosek
from   mosek.fusion import *
from 	 mosek import *
import numpy as np
import time

startTime = time.time()

eps = 1E-10

with Model("Plain Flag Algebra") as M:
        
	# Setting up the variables
	x = M.variable("x", Domain.greaterThan(0.))
	
	print("Reading in A.")
	with open("plainFlagAlgebra1.txt", "r") as file1:
		line = file1.readline()
		dimension = [int(i) for i in line.split(" ") if i.strip()]
		
		firstIndex = 0
		Matrices = [[[0]*dimension[2] for k in range(dimension[2])] for j in range(dimension[1])]
		
		A = []
		for i in range(dimension[0]):
			A.append(M.variable("A"+str(i), Domain.inPSDCone(dimension[i+2]))) 
			
		constr = []
		
		for line in file1.readlines():
			temp = [float(i) for i in line.split(" ") if i.strip()]
			
			if firstIndex != temp[0]:
				print("Index is ", int(temp[0]), " out of ", dimension[0])
				if firstIndex == 0:
					for i in range(dimension[1]):
						tempMat = Matrix.dense(Matrices[i])
						constr.append(Expr.add(Expr.dot(tempMat, A[0]), x))
						
				else:
					for i in range(dimension[1]):
						tempMat = Matrix.dense(Matrices[i])
						constr[i] = Expr.add(constr[i],Expr.dot(tempMat, A[firstIndex]))
			
				Matrices = [[[0]*dimension[int(temp[0])+2] for k in range(dimension[int(temp[0])+2])] for j in range(dimension[1])]
				firstIndex = int(temp[0])
			
			Matrices[int(temp[1])][int(temp[2])][int(temp[3])] = temp[4]-eps
			Matrices[int(temp[1])][int(temp[3])][int(temp[2])] = temp[4]-eps
			
	#Add final constraints
	for i in range(dimension[1]):
		temp = Matrix.dense(Matrices[i])
		constr[i] = Expr.add(constr[i],Expr.dot(temp, A[firstIndex]))
	
	print("Reading in B.")		
	with open("plainFlagAlgebra2.txt", "r") as file1:
		for line in file1.readlines():
			B = [float(i)-eps for i in line.split(" ") if i.strip()]
	
	print("Reading in C.")
	C = []
	numKnown = 0
	with open("plainFlagAlgebra3.txt", "r") as file1:
		for line in file1.readlines():
			test = [float(i) for i in line.split(" ") if i.strip()]
			C.append(test)
			numKnown = numKnown + 1

	if dimension[1] != len(B) :
		print("Something wrong with .txt files.")
	

	#Creating some variable to bound dual solutions
	#To see this look at how to bound duals in LP
	y = []
	
	for i in range(numKnown):
		y.append(M.variable("y"+str(i), Domain.greaterThan(0.)))
	
	#Add constraint on duals to bound subgraph densities
	forObj = Expr.mul(-C[0][0]-eps,y[0])
	for i in range(numKnown) :
		if i != 0 :
			forObj = Expr.add(forObj,Expr.mul(-C[i][0]-eps,y[i]))
			
		for j in range(dimension[1]) : 
			constr[j] = Expr.add(constr[j],Expr.mul(y[i],-C[i][j+1]+eps))
		
	#Add contraint that all duals sum <= 1
	#z = M.variable("z2", Domain.greaterThan(0.))
	
	#forObj = Expr.add(forObj,Expr.neg(z))
	
	#for i in range(dimension[1]) :
		#constr[i] = Expr.add(constr[i],Expr.neg(z))
		
	X = []

	for i in range(dimension[1]) :
		#c = M.constraint("c"+str(i),constr[i], Domain.lessThan(B[i])) #Minimize
		c = M.constraint("c"+str(i),constr[i], Domain.lessThan(1.-B[i])) #Maximize
		X.append(c)	
			
	# Objective
	# Bounding duals creates new variables, which are in the objective function
	M.objective(ObjectiveSense.Maximize, Expr.add(x,forObj))
	
	M.setLogHandler(sys.stdout)
	M.writeTask("sdo3.ptf")
	M.writeTask("dump.task.gz")
	M.writeTask("data.opf")
	M.solve()

	index = 1
	
	print("Printing all non-zero (> 1E-5) duals:")
	for c in X:
		if(c.dual() > 0.00001): 
			print("For graph ", index, " we have a dual value of: ", c.dual())
		index = index + 1
	#Minimize
	#print("Solution is: ", M.dualObjValue())
	
	#Maximize
	print("Solution is: ", 1-M.dualObjValue())

print("Total running time: ", time.time()-startTime, " seconds. ")
