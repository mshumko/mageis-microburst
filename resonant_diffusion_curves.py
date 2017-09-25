import numpy as np

# Define physical constants
c = 3.0E8 # m/s
mu_0 = 1.26E-6 # H·m−1 or N·A−2
eps0 = 8.85E-12 # F/m
m_p = 1.7E-27 # kg
m_e = 9.1E-31 # kg
q_e = -1.6E-19 # C

B0 = 31.2E-6 # Earth's magnetic field strength in Tesla from Schultz and Lanzerotti.

# MagB is from Eq. 1.23 in Schultz and Lanzerotti for a dipole field.
magB = lambda mlat, L: (B0/L**3)*np.sqrt(1 + 3*np.power(np.sin(np.deg2rad(
    mlat)), 2))/np.cos(np.deg2rad(mlat))**6

# Cyclotron frrequency
wce = lambda λ, Ldip: np.abs(q_e)*magB(λ, Ldip)/m_e 
# Electron number density. Currently constant, but can assume a complex function.
n_e = lambda n0, mlat = None: n0 
# Plasma frequency
wpe = lambda n0, mlat = None: np.sqrt(n_e(n0, mlat)*q_e**2/(m_e*eps0))
# Chorus |k| in a cold plama 
magk = lambda w, n0, λ, Ldip: (w/c)*np.sqrt(1 - wpe(n0, λ)**2/(
    w*(w - wce(λ, Ldip))))

def p(vPerp, vParallel):
    """
    NAME:    p(vPerp, vParallel)
    USE:     Map velocity to momentum space with relativstic effects accounted for.
    INPUT:   Parallel and perpendicular velocity in m/s.
    AUTHOR:  Mykhaylo Shumko
    RETURNS: perpendicular and parallel momenta in normalized momnetum 
             (Momentum in SI units)/(me*c)
    MOD:     2017-09-25
    """
    validInd = np.where(np.isfinite(vPerp))
    v = np.sqrt(np.power(vPerp[validInd], 2) + np.power(vParallel[validInd], 2))
    g = 1/np.sqrt(1 - v**2/c**2)
    pPerp = vPerp[validInd]*g/c
    pParallel = vParallel[validInd]*g/c
    return pPerp, pParallel

def resCurveVperp(vParallel, w, n0, mlat, L, n = 1):
    """
    NAME:    resCurveVperp(vParallel, w, n0, mlat, L, n = 1)
    USE:     Uses the resonance condition to return the perpendicular velocity.
             NOTE: For electron-chorus scattering, vParallel must be negative!
    INPUT:   REQUIRED:
                vParallel - Parallel velocity in m/s
                w - Wave frequency in Rad/s
                n0 - Electron number density in #/m^3
                mlat - Magnetic latitude (for calculating |B_dipole|)
                L - Dipole L shell (for calculating |B_dipole|)
             OPTIONAL:
                n = 1 - Resonance harmonic (default to cyclotron resonance)
    AUTHOR:  Mykhaylo Shumko
    RETURNS: Perpendicular velocity in m/s.
    MOD:     2017-09-25
    """
    A = (c*(w - vParallel*magk(w, n0, mlat, L))/(n*wce(mlat, L)))**2
    return np.sqrt(c**2 - vParallel**2 - A)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    print('Recreating resonance curves from Figure 6 in Meredith et al 2002.')
    ### CITATION ###
    # N. P. Meredith, R. B. Horne, D. Summers, R. M. Thorne, R. H. A. 
    # Iles, et al.. Evidence for acceleration of outer zone electrons to 
    # relativistic energies by whistler mode chorus. Annales
    # Geophysicae, European Geosciences Union, 2002, 20 (7), pp.967-979.

    vParallel = c*np.linspace(0, -0.99, num = 100) # m/s
    mlat = 0
    L = 5.7 #4
    n0 = 0.5E6 #60.5E6 # Density (# m^-3)

    # Calculate and plot resutls
    vPerp = resCurveVperp(vParallel, 0.1*wce(mlat, L), n0, mlat, L)
    pPerp, pParallel = p(vPerp, vParallel)
    plt.plot(pPerp, -pParallel, label = r'$0.1 \ \Omega_{ce}$')

    vPerp = resCurveVperp(vParallel, 0.3*wce(mlat, L), n0, mlat, L)
    pPerp, pParallel = p(vPerp, vParallel)
    plt.plot(pPerp, -pParallel, label = r'$0.3 \ \Omega_{ce}$')

    vPerp = resCurveVperp(vParallel, 0.6*wce(mlat, L), n0, mlat, L)
    pPerp, pParallel = p(vPerp, vParallel)
    plt.plot(pPerp, -pParallel, label = r'$0.6 \ \Omega_{ce}$')

    plt.ylim((-1, 4))
    plt.xlim((0, 4))
    plt.xlabel(r'$p_{\perp}/m_e c$')
    plt.ylabel(r'$p_{||}/m_e c$')
    plt.title(r'Resonance Curves for $\omega_{pe} / \omega_{ce} = $' + str(round(wpe(n0, mlat)/wce(mlat, L), 2)))
    plt.legend()
    plt.show()

    vParallel_res = c*np.linspace(0, -0.99, num = 1000)
    mlat = 0
    L = 5.7
    n0 = 0.5E6 # Density at the time

    vPerp_res = resCurveVperp(
        vParallel_res, 0.1*wce(mlat, L), n0, mlat, L)
    pPerp_res, pParallel_res = p(vPerp_res, vParallel_res)
    print(pPerp_res, pParallel_res)
    plt.plot(pPerp_res, -pParallel_res, 'k--', label=r'$0.9 \ \Omega_{ce}$')
    plt.show()
