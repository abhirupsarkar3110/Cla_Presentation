import numpy as np
from nrminv import _setup, _factor

def nrminv_final1(a, b, c):
    INF = np.inf
    n, A, B, C = _setup(a, b, c)
    with np.errstate(all="ignore"):
        f = _factor(n, A, B, C)
        if f is None:
            return INF
        Dp, Dm, Lp, Um = f
        if Dp[n] == 0 or Dm[1] == 0:
            return INF
        for i in range(2, n):
            if Dp[i-1] == 0 and Dm[i+1] == 0:
                return INF
        
        delta = np.zeros(n+3)
        su = np.zeros(n+3)
        sl = np.zeros(n+3)
        
        delta[1] = 1.0 / Dm[1]
        for i in range(1, n):
            if Dp[i] == 0 or 1.0 / Dm[i+1] == 0: # pivot zero / infinite
                delta[i+1] = 0.0
            elif Dm[i+1] == 0: # eq. (6.2)
                delta[i+1] = -Dp[i] / (B[i] * C[i])
            elif 1.0 / Dp[i] == 0: # eq. (6.1)
                delta[i+1] = 1.0 / Dm[i+1]
            elif delta[i] == 0: # eq. (8.1), underflow
                delta[i+1] = 1.0 / (Dm[i+1] - B[i] * C[i] / Dp[i])
            else: # eq. (4.11)
                delta[i+1] = delta[i] / Dm[i+1] * Dp[i]
                
            if 1.0 / delta[i+1] == 0:
                return INF
                
        for i in range(1, n): # upper column sums
            if su[i] + abs(delta[i]) == 0: # eq. (6.3)
                su[i+1] = (su[i-1] + abs(delta[i-1])) * abs(C[i-1] / B[i]) + abs(1.0 / B[i])
            else: # eq. (4.12)
                su[i+1] = (su[i] + abs(delta[i])) * abs(Um[i])
                
            if 1.0 / su[i+1] == 0:
                return INF
                
        for i in range(n, 1, -1): # lower column sums
            if sl[i] + abs(delta[i]) == 0: # eq. (6.4)
                sl[i-1] = (sl[i+1] + abs(delta[i+1])) * abs(B[i] / C[i-1]) + abs(1.0 / C[i-1])
            else: # eq. (4.13)
                sl[i-1] = (sl[i] + abs(delta[i])) * abs(Lp[i-1])
                
            if 1.0 / sl[i-1] == 0:
                return INF
                
        gamma = 0.0
        for i in range(1, n+1):
            if su[i] + sl[i] + abs(delta[i]) > gamma:
                gamma = su[i] + sl[i] + abs(delta[i])
                
        return gamma




#To check if code is running properly, you can uncomment the following test code and run the script. It will calculate the norm for a sample tridiagonal matrix defined by the vectors a, b, and c.

# if __name__ == "__main__":
#     a = [2.0, 2.0, 2.0, 2.0, 2.0]  # Main diagonal (length n)
#     b = [-1.0, -1.0, -1.0, -1.0]   # Sub-diagonal (length n-1)
#     c = [-1.0, -1.0, -1.0, -1.0]   # Super-diagonal (length n-1)
    
#     print("Running nrminv_final1 test...")
#     try:
#         result = nrminv_final1(a, b, c)
#         print(f"Success! Calculated norm: {result}")
#     except Exception as e:
#         print(f"Error occurred: {e}")