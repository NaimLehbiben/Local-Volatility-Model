import numpy as np
from scipy.optimize import minimize


def sabr_vol(T, K, F0, sigma0, alpha, rho, beta=0.3) :
    """
    Hagan's 2002 SABR log-normal vol expansion.
    The strike K can be a scalar or an array, the function will return an array
    of log-normal vols.
    """
    eps = 1e-04
    fk_beta = (F0 * K)**(1 - beta)
    log_fk = np.log(F0/K)

    a = ((1 - beta)**2 / 24) * (sigma0**2 / fk_beta)
    b = 0.25 * rho * beta * alpha * sigma0 / fk_beta**0.5
    c = (2 - 3 * rho**2) * alpha**2 / 24
    d = fk_beta**0.5
    e = (1 - beta)**2 * log_fk**2 / 24
    f = (1 - beta)**4 * log_fk**4 / 1920
    z = alpha * fk_beta**0.5 * log_fk / sigma0

    # For close points to ATM, z and x(z) become close to zero and should be removed for the Hagan's formula
    num = sigma0 * (1 + (a + b + c) * T)
    num = np.where(abs(z) > eps, num * z, num)
    den = d * (1 + e + f)
    den = np.where(abs(z) > eps, den * x(rho, z), den)

    return num / den

def x(rho, z) :
    """
    Return function x used in Hagan's 2002 SABR log-normal vol expansion.
    """
    num = (1 - 2 * rho * z + z**2)**0.5 + z - rho
    den = 1 - rho
    return  np.log (num / den)


def sabr_calibration(T, F0, strikes, vol_market, beta = 0.3) :
    """
    SABR model calibration : Sigma0, Alpha, and Rho.
    """
    parameters0 = np.array([vol_market[len(vol_market) // 2], 0.2, -0.2]) # vol_market[len(vol_market) // 2] is the closest point vol to ATM, 0.2 and -0.2 were chosen arbitrary
    res = minimize(obj_func, x0=parameters0, args=(T, F0, strikes, vol_market, beta), bounds=((0.01, None), (0.01, None), (-0.99, 0.99)))
    return res.x

def obj_func(parameters, T, F0, strikes, vol_market, beta) :
    """
    The Sum Squared Error function : Used as the objective function in the SABR Calibration.
    """
    sabr_volatilities = sabr_vol(T, strikes, F0, parameters[0], parameters[1], parameters[2], beta)
    sse = 0
    for i in range(len(strikes)) :
        sse += (sabr_volatilities[i] - vol_market[i])**2
    return sse

def implied_vol_sabr(T, K) :
    """
    SABR Implied vol : Time dependent Monotonic Cubic interpolation of the SABR parameters.
    """
    F0 = S0 * np.exp(zc_rates_interp(T) * T)
    sigma0 = sigma_interp(T)
    alpha = alpha_interp(T)
    rho = rho_interp(T)
    return sabr_vol(T, K, F0, sigma0, alpha, rho)

def local_vol(t, S) :
    """
    Local Vol from Dupire's formula with no Dividends.
    """
    t = np.where(t == 0, 1e-04, t) # to avoid division by 0
    eps_t = 1e-05
    eps_S = eps_t * S
    vol = implied_vol_sabr(t, S)
    sqrt_t = np.sqrt(t)
    zc_rate_t = zc_rates_interp(t)
    r_t = zc_rate_t + t * (zc_rates_interp(t + eps_t) - zc_rate_t) / eps_t # r(T) = zc_r(T) + T * (dzc_r(T) / dT)
    dvol_dT = (implied_vol_sabr(t + eps_t, S) - vol) / eps_t
    dvol_dK = (implied_vol_sabr(t, S + eps_S) - vol) / eps_S
    d2vol_dK2 = (implied_vol_sabr(t, S + eps_S) + implied_vol_sabr(t, S - eps_S) - 2 * vol ) / eps_S**2
    d1 = (np.log(S0 / S) + zc_rate_t * t + 0.5 * vol**2 * t) / (vol * sqrt_t)
    num = 2 * dvol_dT + vol / t + 2 * S * r_t * dvol_dK
    den = S**2 * (d2vol_dK2 - d1*sqrt_t*dvol_dK**2 + (1/(S*sqrt_t) + d1*dvol_dK)**2 / vol)
    return np.sqrt(np.where(num/den >= 0, num/den, 0))