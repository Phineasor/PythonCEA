# Bisect Function for root finding zeros of a defined function in python, the function needs to be only dependent on x
def Bisect(f, x1, x2, tol) -> float:
    Error = 1
    xm = x2
    i = 0
    while (Error > tol) & (i < 500):
        i += 1
        xmP = xm
        xm = (x1 + x2) / 2

        if f(xm) * f(x1) < 0:
            x2 = xm
            Error = abs(xmP - xm) / xm
        else:
            x1 = xm
            Error = abs(xmP - xm) / xm
    return xm
