import numpy as np
from nrminv import _setup, _factor

INF = np.inf
def nrminv_final2(a, b, c):
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
        
        for i in range(1, n): # eq. (5.6)
            delta[i+1] = 1.0 / (Dm[i+1] - B[i] * C[i] / Dp[i])
            if 1.0 / delta[i+1] == 0:
                return INF
                
        for i in range(1, n): # upper column sums
            DET = Dm[i+1] * A[i] - B[i] * C[i]
            if 1.0 / Dm[i+1] == 0 or abs(Dm[i+1] * A[i]) >= abs(DET):
                su[i+1] = (su[i] + abs(delta[i])) * abs(Um[i]) # ordinary formula
            else: # block pivot, eq. (8.5)
                su[i+1] = ((su[i-1] + abs(delta[i-1])) * abs(C[i-1] * C[i] / DET) 
                           + abs(C[i] / (Dm[i+1] * Dp[i] - B[i] * C[i])))
            if 1.0 / su[i+1] == 0:
                return INF
                
        for i in range(n, 1, -1): # lower column sums
            DET = Dp[i-1] * A[i] - B[i-1] * C[i-1]
            if 1.0 / Dp[i-1] == 0 or abs(Dp[i-1] * A[i]) >= abs(DET):
                sl[i-1] = (sl[i] + abs(delta[i])) * abs(Lp[i-1])
            else:
                sl[i-1] = ((sl[i+1] + abs(delta[i+1])) * abs(B[i-1] * B[i] / DET) 
                           + abs(B[i-1] / (Dp[i-1] * Dm[i] - B[i-1] * C[i-1])))
            if 1.0 / sl[i-1] == 0:
                return INF
                
        gamma = 0.0
        for i in range(1, n+1):
            if su[i] + sl[i] + abs(delta[i]) > gamma:
                gamma = su[i] + sl[i] + abs(delta[i])
                
        return gamma




#To check if code is running properly, you can uncomment the following test code and run the script. It will calculate the norm for a sample tridiagonal matrix defined by the vectors a, b, and c.


# if __name__ == "__main__":
#     a = [2, 2, 2, 2, 2]    # Main diagonal
#     b = [-1, -1, -1, -1]   # Sub-diagonal
#     c = [-1, -1, -1, -1]   # Super-diagonal
    
#     print("Running nrminv_final2 test...")
#     try:
#         result = nrminv_final2(a, b, c)
#         print(f"Success! Calculated norm: {result}")
#     except Exception as e:
#         print(f"Error occurred: {e}")