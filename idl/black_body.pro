function b_planck, T, nu
  COMMON CONSTANTS
  fact = h_planck*nu/(k_B * T)
  bb = 2.0/(speed_of_light^2) * h_planck * nu*nu*nu / (exp(fact)-1.0) 
  return, bb
end

r_star = 10.0 * r_sol
d = 10.68*pc
T = 40000.0
nu = 3.288e15
E_ion = 13.6 * erg_eV
nu = E_ion / h_planck
print, "frequency: ",nu
dilu = (r_star/d)^2
flux_bol   = dilu * sigma_SB * T^4
planck_bol = flux_bol / !pi
i_nu = b_planck(40000.0,nu)
nr_of_photons = i_nu / (h_planck*nu)
print, "nr of photons:", nr_of_photons

end

