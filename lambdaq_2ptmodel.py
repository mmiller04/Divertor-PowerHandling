### file used to calculate 2pt model
import numpy as np


def main():

	# these are the parameters needed for calculation - will change over time as design progresses

	R0_m = 3.3
	a_m = 1.1
	eps = a_m/R0_m

	# engineering parameters
	Psol_MW = 25
	Bt_T = 9.2
	q95 = 4.755 # 7.2
	Ip_MA = 7.8
	fGW = 0.67 # ARCH was 0.7375

	mu0 =  4e-7*np.pi
	Bp_T = 1.5 #mu0*(Ip_MA*1e6)/(2*R0_m)

	# plasma parameters
	ne0_m3 = 4.3
	Te0_keV = 20
	pe0_MPa = 1.4

	# estimates for quick and dirty calculations of volume-averaged plasma quantities
	ne_coef = 3/4
	Te_coef = 1/2

	nevol_m3 = ne_coef*ne0_m3
	Tevol_keV = Te_coef*Te0_keV
	pevol_atm = ne_coef*Te_coef*(pe0_MPa * 9.86923)

	pvol_atm = 2*pevol_atm # assuming pe = pi - may need to be checked


	# calcluate lambda qs!

	lq_eich14 = lambdaq_eich14(Bp_T)
	lq_eich15 = lambdaq_eich15(Psol_MW, R0_m, Bp_T, eps)
	lq_brunner = lambdaq_brunner(pvol_atm)
	lq_brunnerL = lambdaq_brunner(Bp_T, mode='L')
	lq_brunnerI = lambdaq_brunner(Bp_T, mode='I')
	lq_brunnerH = lambdaq_brunner(Bp_T, mode='H')
	lq_horacek = lambdaq_horacek(Bp_T, q95, fGW)



	print('Eich14: {} mm'.format(lq_eich14))
	print('Eich15: {} mm'.format(lq_eich15))
	print('Brunner: {} mm'.format(lq_brunner))
	print('Brunner-Eich L-mode: {} mm'.format(lq_brunnerL))
	print('Brunner-Eich I-mode: {} mm'.format(lq_brunnerI))
	print('Brunner-Eich H-mode: {} mm'.format(lq_brunnerH))
	print('Horacek: {} mm'.format(lq_horacek))



def two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, model='eich'):


	R_sep_m = R0_m + a_m
	eps = a_m/R0_m

	# coefficients for heat conduction by electrons or H ions

	k0_e = 2000
	k0_i = 2000
	gamma = 7 # sheath heat flux transmission coefficient


	return Te_sep



def lambdaq_eich14(Bp_T): # scaling #14 from Eich et al.
	
	C0 = 0.63
	a = -1.19

	lambdaq_mm = C0*(Bp_T**a)

	return lambdaq_mm



def lambdaq_eich15(Psol_MW, R0_m, Bp_T, eps): # scaling #15 from Eich et al.  
	
	C0 = 1.35
	a = -0.02
	b = 0.04
	c = -0.92
	d = 0.42

	lambdaq_mm = C0*(Psol_MW**a)*(R0_m**b)*(Bp_T**c)*(eps**d)

	return lambdaq_mm
	


def lambdaq_brunner(pvol_atm): # scaling from Brunner et al.

	C0 = 0.91
	a = -0.48

	lambdaq_mm = C0*(pvol_atm**a)

	return lambdaq_mm


def lambdq_brunner_highfield(Bp_T, mode='L'): # scaling from Brunner et al. for C-Mod across confinement regimes

	if mode == 'L':
		C0 = 1.37
		a = -0.74
	if mode == 'L':
		C0 = 0.95
		a = -0.57
	if mode == 'L':
		C0 = 0.76
		a = -0.96

	lambdaq_mm = C0*(Bp_T**a)

	return lambdaq_mm


def lambdaq_horacek(Bp_T, q95, fGW): # l-mode scaling from Horacek et al.

	C0 = 8.39
	a = -0.36
	b = 0.55
	c = 0.92

	lambdaq_mm = C0*(Bp_T**a)*(q95**b)*(fGW**c)

	return lambdaq_mm


if __name__ == '__main__':

	main()






