### file used to calculate 2pt model
import numpy as np


def main():

	# these are the parameters needed for calculation - will change over time as design progresses

	R0_m = 4
	a_m = 1
	eps = a_m/R0_m

	# engineering parameters
	Psol_MW = 25
	Bt_T = 11.025
	q95 = 4.755 # 7.2
	# p95 = 
	Ip_MA = 9.74
	fGW = 0.67 # ARCH was 0.7375

	mu0 =  4e-7*np.pi
	Bp_T = 1.5 #mu0*(Ip_MA*1e6)/(2*R0_m)

	# plasma parameters
	Te0_keV = 21.5
	pe0_MPa = 1.8

	ne0_m3 = pe0_MPa*1e6 / (Te0_keV*1e3 / 6.242e18)

	# estimates for quick and dirty calculations of volume-averaged plasma quantities
	ne_coef = 3/4
	Te_coef = 1/2

	nevol_m3 = ne_coef*ne0_m3
	Tevol_keV = Te_coef*Te0_keV
	pevol_atm = ne_coef*Te_coef*(pe0_MPa * 9.86923)

	pvol_atm = 2*pevol_atm # assuming pe = pi - may need to be checked


	# calcluate lambda qs!

	# lq_eich14 = lambdaq_eich14(Bp_T)
	# lq_eich15 = lambdaq_eich15(Psol_MW, R0_m, Bp_T, eps)
	lq_brunner = lambdaq_brunner(pvol_atm)
	lq_brunnerL = lambdaq_brunner_highfield(Bp_T, mode='L')
	lq_brunnerI = lambdaq_brunner_highfield(Bp_T, mode='I')
	lq_brunnerH = lambdaq_brunner_highfield(Bp_T, mode='H')
	lq_horacek = lambdaq_horacek(Bp_T, q95, fGW)
	lq_scarabosio = lambdaq_scarabosio(Bt_T, Psol_MW, q95, R0_m)
	# lq_silvagni = lambdaq_silvagni(p95)

	Tesep_brunner = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lq_brunner)
	Tesep_brunnerL = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lq_brunnerL)
	Tesep_brunnerI = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lq_brunnerI)
	Tesep_brunnerH = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lq_brunnerH)
	Tesep_scarabosio = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lq_scarabosio)
	# Tesep_silvagni = two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lambdaq_eich14)


	# print('Eich14: {} mm, {} eV'.format(lq_eich14, Tesep_eich14))
	# print('Eich15: {} mm, {} eV'.format(lq_eich15, Tesep_eich15))
	print('Brunner: {} mm, {} eV'.format(lq_brunner, Tesep_brunner))
	print('Brunner-Eich L-mode: {} mm, {} eV'.format(lq_brunnerL, Tesep_brunnerL))
	print('Brunner-Eich I-mode: {} mm, {} eV'.format(lq_brunnerI, Tesep_brunnerI))
	print('Brunner-Eich H-mode: {} mm, {} eV'.format(lq_brunnerH, Tesep_brunnerH))
	# print('Horacek: {} mm, {} eV'.format(lq_horacek, Tesep_horacek))
	print('Scarabioso: {} mm, {} eV'.format(lq_scarabosio, Tesep_scarabosio))


def two_point_model(R0_m, a_m, Psol_MW, Bp_T, Bt_T, q95, lambdaq_mm):


	R_sep_m = R0_m + a_m

	# coefficients for heat conduction by electrons or H ions

	k0_e = 2000
	k0_i = 2000
	gamma = 7 # sheath heat flux transmission coefficient

	L_par = np.pi*R_sep_m*q95

	q_par_MW_m2 = 1/2*Psol_MW / (2*np.pi*R_sep_m*(lambdaq_mm * 1e-3)) * np.hypot(Bt_T, Bp_T) / Bp_T


	Tu_eV = ((7/2) * (q_par_MW_m2*1e6) * L_par / (2*k0_e))**(2/7)

	return Tu_eV



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


def lambdaq_brunner_highfield(Bp_T, mode='L'): # scaling from Brunner et al. for C-Mod across confinement regimes

	if mode == 'L':
		C0 = 1.37
		a = -0.74
	if mode == 'I':
		C0 = 0.95
		a = -0.57
	if mode == 'H':
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


def lambdaq_scarabosio(Bt_T, Psol_MW, q95, R0_m):

	C0 = 1.58
	a = -0.4
	b = 0.13
	c = 0.73
	d = 0.26

	lambdaq_mm = C0*(Bt_T**a)*(Psol_MW**b)*(q95**c)*(R0_m**0.26)

	return lambdaq_mm


def lambdaq_silvagni(p95): # l-mode scaling from silvagni

	C0 = 2.45
	a = -0.34

	lambdaq_mm = C0*(p95**a)

	return lambdaq_mm


if __name__ == '__main__':

	main()






