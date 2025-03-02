from scipy.stats import norm
import math

def black_scholes(flavor, K, T, S, r, vol, q=0) :
    """
    The Black-Scholes price.
    """
    v2T = vol**2 * T
    d1 = (math.log(S/K) + (r - q) * T + v2T / 2) / v2T**0.5
    d2 = d1 - v2T**0.5
    phi = 1 if flavor.upper() == "CALL" else  -1
    return phi * (S * math.exp(-q * T) * norm.cdf(phi * d1) - K * math.exp(-r * T) * norm.cdf(phi * d2))

def vega_bs(K, T, S, r, vol, q=0) :
    """
    The Black-Scholes Vega greek : The derivative of the option value with respect to the volatility of the underlying asset.
    """
    v2T = vol**2 * T
    d2 = (math.log(S/K) + (r - q) * T - v2T / 2) / v2T**0.5
    return K * math.exp(-r * T) * norm.pdf(d2) * T**0.5

import sys

def newton_raphson(market_price, init_vol, flavor, K, T, S, r, q=0):
    """
    The Newton-Raphson algorithm : Implied volatilities from market prices
    """
    eps = 1e-07
    vol = init_vol
    func = black_scholes(flavor, K, T, S, r, vol, q) - market_price
    func_deriv = vega_bs(K, T, S, r, vol, q)
    error = sys.exit("The algorithm failed to converge. Please review the input data.") if func_deriv < 0.0001 else "No error."
    next_vol = vol - func / func_deriv
    error = sys.exit("The algorithm failed to converge. Please review the input data.") if next_vol < 0.0001 else "No error."
    nb_iteration = 1
    while abs(next_vol - vol) > eps :
        vol = next_vol
        func = black_scholes(flavor, K, T, S, r, vol, q) - market_price
        func_deriv = vega_bs(K, T, S, r, vol, q)
        error = sys.exit("The algorithm failed to converge. Please review the input data.") if func_deriv < 0.0001 else "No error."
        next_vol = vol - func / func_deriv
        error = sys.exit("The algorithm failed to converge. Please review the input data.") if next_vol < 0.0001 else "No error."
        nb_iteration += 1
        error = sys.exit("The algorithm failed to converge. Please review the input data.") if nb_iteration > 300 else "No error."

    return next_vol