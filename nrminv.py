import numpy as np

def _setup(a, b, c):
    n = len(a)
    A = np.zeros(n+3)
    B = np.zeros(n+3)
    C = np.zeros(n+3)
    
    A[1:n+1] = a
    B[1:n] = b
    C[1:n] = c
    return n, A, B, C

def _factor(n, A, B, C, check=True):
    Dp = np.zeros(n+3)
    Dm = np.zeros(n+3)
    Lp = np.zeros(n+3)
    Um = np.zeros(n+3)
    
    Dp[1] = A[1]
    for i in range(1, n):
        if check and Dp[i] == 0 and B[i] * C[i] == 0:
            return None
        Lp[i] = B[i] / Dp[i]
        Dp[i+1] = A[i+1] - C[i] * Lp[i]
        
    Dm[n] = A[n]
    for i in range(n-1, 0, -1):
        if check and Dm[i+1] == 0 and B[i] * C[i] == 0:
            return None
        Um[i] = C[i] / Dm[i+1]
        Dm[i] = A[i] - B[i] * Um[i]
        
    return Dp, Dm, Lp, Um

def nrminv(a, b, c):
    n, A, B, C = _setup(a, b, c)
    with np.errstate(all="ignore"):
        Dp, Dm, Lp, Um = _factor(n, A, B, C, check=False)
        delta = np.zeros(n+3)
        su = np.zeros(n+3)
        sl = np.zeros(n+3)
        
        delta[1] = 1.0 / Dm[1]
        
        for i in range(1, n): # eq. (4.11)
            delta[i+1] = delta[i] * Dp[i] / Dm[i+1]
            
        for i in range(n-1, 0, -1): # eq. (4.12)
            su[i] = (su[i+1] + abs(delta[i+1])) * abs(Um[i])
            
        for i in range(2, n+1): # eq. (4.13)
            sl[i] = (sl[i-1] + abs(delta[i-1])) * abs(Lp[i-1])
            
        return max(su[i] + sl[i] + abs(delta[i]) for i in range(1, n+1))