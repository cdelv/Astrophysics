import numpy as np
from scipy.integrate import odeint
import math
import matplotlib.pyplot as plt

import sympy
from einsteinpy.symbolic import MetricTensor, RiemannCurvatureTensor, predefined

coord = sympy.symbols('r theta phi t')

def main():

	##### black whole mass
	M = 1000

	#Initial conditions   r, phi.
	z0=[50,0]
	tmax = 20
	steps = 100

	#Mass of the particle
	m = 1 

	#Energy of the particle
	E = 10

	#Orbital angular momentum
	j = 1

	#Spin Angular momentum
	s = 1

	#Create the metric tensor
	Metric_Tensor=Create_Metric_Tensor(M)
	up_Metric_Tensor = Metric_Tensor.change_config('uu') #metric tensor with 2 up indices

	print(Metric_Tensor)
	
	#Create the Riemann curvature tensor from the metric
	Riemann_Tensor = Create_Riemann_Tensor(Metric_Tensor)

	## We are always in the plain. Evaluate theta = pi/2 to speed up the code
	## We are not intrested on the complete tensor, just some slices
	R_30 = Riemann_Tensor[3][0].subs(coord[1], math.pi/2)
	R_31 = Riemann_Tensor[3][1].subs(coord[1], math.pi/2)
	R_10 = Riemann_Tensor[1][0].subs(coord[1], math.pi/2)
	R_13 = Riemann_Tensor[1][3].subs(coord[1], math.pi/2)
	Metric_Tensor = Metric_Tensor.subs(coord[1], math.pi/2)
	up_Metric_Tensor = up_Metric_Tensor.subs(coord[1], math.pi/2)

	## Calculate analitic expression of gama and its derivative
	gamma, Part_gamma = calc_gamma(Metric_Tensor)

	### Calculate analitic Metric partial derivatives
	Part_g00 = sympy.diff(Metric_Tensor[0][0], coord[0])
	Part_g33 = sympy.diff(Metric_Tensor[3][3], coord[0])

	#Total Angular momentum
	J = j + s

	#time span
	t=np.linspace(0,tmax,steps) #start, finish, n of points

	#Solve ODE
	z = odeint(dif_equation, z0, t, args=(Metric_Tensor, up_Metric_Tensor, 
		R_30, R_31, R_10, R_13, m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33))

	#Plot trayectory
	plot_trayectory(z,t) #r vs phi

def Create_Metric_Tensor(M=10):
	#### Define problem metric #########
	#### schwarzschild metric ##########
	#list2d = [[0 for i in range(4)] for i in range(4)]
	#list2d[0][0] = -(1-2*M/coord[0])
	#list2d[1][1] = 1 / (1 - 2*M / coord[0])
	#list2d[2][2] = coord[0]**2
	#list2d[3][3] = coord[0]**2*sympy.sin(coord[1])**2
	#Metric = MetricTensor(list2d, coord)
	#Metric.tensor()

	########### premade schwarzschild metric #################### c, G, a, M
	Metric = predefined.janis_newman_winicour.JanisNewmanWinicour(1, 1, 1, M)
	return Metric  ### Metric_{ab}  a, b in (0 to 4), 2 indices down

def Create_Riemann_Tensor(Metric):
	#### Calculate Riemann Tensor from the metric tensor #########
	Riemann = RiemannCurvatureTensor.from_metric(Metric)
	Riemann.tensor()
	return Riemann

def calc_gamma(Metric):

	#calculate analitic expresion of gamma and its derivative
	gamma = sympy.sqrt(-Metric[1][1]*Metric[2][2]*Metric[3][3])
	Part_gamma = sympy.diff(gamma, coord[0])

	return gamma, Part_gamma

def calc_momentum(up_Metric, gamma, Part_g00, Part_g33, m, E, s, J):

	### Convert tu number components of metric tensor
	g_up_00 = sympy.N(up_Metric[0][0])
	g_up_11 = sympy.N(up_Metric[1][1])
	g_up_33 = sympy.N(up_Metric[3][3])

	P0 = -(2*m*gamma*(2*m*gamma*E+s*J*Part_g00))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	P3 = (2*m*gamma*(2*m*gamma*J+s*E*Part_g33))/(4*m**2*gamma**2+s**2*Part_g00*Part_g33)
	P1_2 = -(m**2+g_up_00*P0**2+g_up_33*P3**2)/(g_up_11)

	return sympy.N(P0), math.sqrt(sympy.N(P1_2)), sympy.N(P3)
	
def calc_spin_tensor(gamma, P, m, s):
	num = -s/(m*gamma)
	S_13 = P[0]*num
	S_01 = P[2]*num
	S_03 = P[1]*num

	return S_13, S_01, S_03

def calc_drdt(R_30, R_31, S, P, s, m, gamma, Part_gamma):

	### Calculate R_30ab * S^ab
	prod1 = R_30[1][3]*S[0] + R_30[0][1]*S[1] + R_30[0][3]*S[2]

	### Calculate R_31ab * S^ab
	prod2 = R_31[1][3]*S[0] + R_31[0][1]*S[1] + R_31[0][3]*S[2]

	fac1 = s/(2*m*gamma)
	fac2 = s/(2*m*gamma**2)

	derivative = (P[1] + fac1*prod1)/(P[0]-P[2]*fac2*Part_gamma-fac1*prod2)

	return derivative

def calc_dpdt(R_10, R_13, S, P, s, m, gamma, Part_gamma, drdt):
	
	### Calculate R_30ab * S^ab
	prod1 = R_10[1][3]*S[0] + R_10[0][1]*S[1] + R_10[0][3]*S[2]

	### Calculate R_31ab * S^ab
	prod2 = R_13[1][3]*S[0] + R_13[0][1]*S[1] + R_13[0][3]*S[2]

	fac1 = s/(2*m*gamma)
	fac2 = s/(2*m*gamma**2)

	derivative = (P[2]-fac2*P[1]*drdt*Part_gamma - fac1*prod1)/(P[0]+fac1*prod2)

	return derivative

def dif_equation(z,t,Metric, up_Metric, R_30, R_31, R_10, R_13, m, E, s, J, gamma, Part_gamma, Part_g00, Part_g33):
	r=z[0]
	phi=z[1]

	### Evaluate the metric tensors and store them in an ndarray
	Metric = Metric.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])
	up_Metric = up_Metric.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])

	### Evaluate Riemann tensor Slices
	R_30 = R_30.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])
	R_31 = R_30.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])
	R_10 = R_30.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])
	R_13 = R_30.subs([(coord[0], r), (coord[2], phi), (coord[3], t)])

	### Evaluate gamma and gamma derivative
	gamma = sympy.N(gamma.subs([(coord[0], r), (coord[2], phi), (coord[3], t)]))
	Part_gamma = sympy.N(Part_gamma.subs([(coord[0], r), (coord[2], phi), (coord[3], t)]))

	### Evaluate metric tensor partial derivatives
	Part_g00 = sympy.N(Part_g00.subs([(coord[0], r), (coord[2], phi), (coord[3], t)]))
	Part_g33 = sympy.N(Part_g33.subs([(coord[0], r), (coord[2], phi), (coord[3], t)]))

	### Calculate momentum
	P = calc_momentum(up_Metric, gamma, Part_g00, Part_g33, m, E, s, J)

	### Calculate non_zero components of the spin tensor
	S = calc_spin_tensor(gamma, P, m, s)

	### Calculate derivatives 
	drdt = calc_drdt(R_30, R_31, S, P, s, m, gamma, Part_gamma)        # r
	dpdt = calc_dpdt(R_30, R_31, S, P, s, m, gamma, Part_gamma, drdt)  #phi

	return drdt, dpdt

def plot_trayectory(z,t):
	r = z[:,0]
	phi = z[:,1]
	x = []
	y = []

	## comvert results to cartesian coordinates for the plot
	for i in range(len(r)):
		x.append(r[i]*math.cos(phi[i]))
		y.append(r[i]*math.sin(phi[i]))
			

	plt.plot(np.asarray(x),np.asarray(y))
	plt.show()

	plt.plot(t,np.asarray(x))
	plt.show()

	plt.plot(t,np.asarray(y))
	plt.show()


if __name__ == '__main__':
	main()