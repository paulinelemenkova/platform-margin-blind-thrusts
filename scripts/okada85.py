import numpy as np

EPS = 1e-14

def _I5(xi, eta, q, dp, nu, R, db):
    X = np.sqrt(xi ** 2 + q ** 2)
    if np.cos(dp) > EPS:
        with np.errstate(divide="ignore", invalid="ignore"):
            I = (1 - 2 * nu) * 2 / np.cos(dp) * np.arctan(
                (eta * (X + q * np.cos(dp)) + X * (R + X) * np.sin(dp)) / (xi * (R + X) * np.cos(dp)))
        I = np.where(xi == 0, 0.0, I)
    else:
        I = -(1 - 2 * nu) * xi * np.sin(dp) / (R + db)
    return I

def _I4(db, eta, q, dp, nu, R):
    if np.cos(dp) > EPS:
        return (1 - 2 * nu) / np.cos(dp) * (np.log(R + db) - np.sin(dp) * np.log(R + eta))
    return -(1 - 2 * nu) * q / (R + db)

def _I3(eta, q, dp, nu, R):
    yb = eta * np.cos(dp) + q * np.sin(dp); db = eta * np.sin(dp) - q * np.cos(dp)
    if np.cos(dp) > EPS:
        return (1 - 2 * nu) * (yb / (np.cos(dp) * (R + db)) - np.log(R + eta)) + \
            np.sin(dp) / np.cos(dp) * _I4(db, eta, q, dp, nu, R)
    return (1 - 2 * nu) / 2 * (eta / (R + db) + yb * q / (R + db) ** 2 - np.log(R + eta))

def _I2(eta, q, dp, nu, R):
    return (1 - 2 * nu) * (-np.log(R + eta)) - _I3(eta, q, dp, nu, R)

def _I1(xi, eta, q, dp, nu, R):
    db = eta * np.sin(dp) - q * np.cos(dp)
    if np.cos(dp) > EPS:
        return (1 - 2 * nu) * (-xi / (np.cos(dp) * (R + db))) - np.sin(dp) / np.cos(dp) * _I5(xi, eta, q, dp, nu, R, db)
    return -(1 - 2 * nu) / 2 * xi * q / (R + db) ** 2

def _atn(xi, eta, q, R):
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.arctan(xi * eta / (q * R))
    return np.where(q != 0, a, 0.0)

def _ux_ss(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    return xi * q / (R * (R + eta)) + _I1(xi, eta, q, dp, nu, R) * np.sin(dp) + _atn(xi, eta, q, R)

def _uy_ss(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    return (eta * np.cos(dp) + q * np.sin(dp)) * q / (R * (R + eta)) + q * np.cos(dp) / (R + eta) + \
        _I2(eta, q, dp, nu, R) * np.sin(dp)

def _uz_ss(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2); db = eta * np.sin(dp) - q * np.cos(dp)
    return db * q / (R * (R + eta)) + q * np.sin(dp) / (R + eta) + _I4(db, eta, q, dp, nu, R) * np.sin(dp)

def _ux_ds(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    return q / R - _I3(eta, q, dp, nu, R) * np.sin(dp) * np.cos(dp)

def _uy_ds(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2)
    return (eta * np.cos(dp) + q * np.sin(dp)) * q / (R * (R + xi)) - \
        _I1(xi, eta, q, dp, nu, R) * np.sin(dp) * np.cos(dp) + np.cos(dp) * _atn(xi, eta, q, R)

def _uz_ds(xi, eta, q, dp, nu):
    R = np.sqrt(xi ** 2 + eta ** 2 + q ** 2); db = eta * np.sin(dp) - q * np.cos(dp)
    return db * q / (R * (R + xi)) - _I5(xi, eta, q, dp, nu, R, db) * np.sin(dp) * np.cos(dp) + \
        np.sin(dp) * _atn(xi, eta, q, R)

def _chin(f, x, p, L, W, q, dp, nu):
    return f(x, p, q, dp, nu) - f(x, p - W, q, dp, nu) - f(x - L, p, q, dp, nu) + f(x - L, p - W, q, dp, nu)

def okada85(e, n, depth, strike, dip, L, W, rake=None, slip=None, ss=None, ds=None, nu=0.25):
    st = np.radians(strike); dp = np.radians(dip)
    if ss is None:
        rk = np.radians(rake); ss = np.cos(rk) * slip; ds = np.sin(rk) * slip
    d = depth + np.sin(dp) * W / 2
    ec = e + np.cos(st) * np.cos(dp) * W / 2
    nc = n - np.sin(st) * np.cos(dp) * W / 2
    x = np.cos(st) * nc + np.sin(st) * ec + L / 2
    y = np.sin(st) * nc - np.cos(st) * ec + np.cos(dp) * W
    p = y * np.cos(dp) + d * np.sin(dp); q = y * np.sin(dp) - d * np.cos(dp)
    ux = -ss / (2 * np.pi) * _chin(_ux_ss, x, p, L, W, q, dp, nu) - ds / (2 * np.pi) * _chin(_ux_ds, x, p, L, W, q, dp, nu)
    uy = -ss / (2 * np.pi) * _chin(_uy_ss, x, p, L, W, q, dp, nu) - ds / (2 * np.pi) * _chin(_uy_ds, x, p, L, W, q, dp, nu)
    uz = -ss / (2 * np.pi) * _chin(_uz_ss, x, p, L, W, q, dp, nu) - ds / (2 * np.pi) * _chin(_uz_ds, x, p, L, W, q, dp, nu)
    return np.sin(st) * ux - np.cos(st) * uy, np.cos(st) * ux + np.sin(st) * uy, uz

def corners(xc, yc, zc, strike, dip, L, W):
    st = np.radians(strike); dp = np.radians(dip)
    us = np.array([np.sin(st), np.cos(st)])
    ud = np.array([np.cos(st), -np.sin(st)]) * np.cos(dp)
    c = np.array([xc, yc])
    tl = c - us * L / 2 - ud * W / 2; tr = c + us * L / 2 - ud * W / 2
    br = c + us * L / 2 + ud * W / 2; bl = c - us * L / 2 + ud * W / 2
    return np.array([tl, tr, br, bl]), zc - np.sin(dp) * W / 2, zc + np.sin(dp) * W / 2
