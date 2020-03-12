from scipy.optimize import newton
import math

# We will use the newton method to find a solution to 
# Y*exp(Y) + 2*N*T*lambda = 0, where Y = -2*(N*lambda + n*nu)*T

def airtime(sf,cr,pl,bw):
    H = 0        # implicit header disabled (H=0) or not (H=1)
    DE = 0       # low data rate optimization enabled (=1) or not (=0)
    Npream = 8   # number of preamble symbol (12.25  from Utz paper)

    if bw == 125 and sf in [11, 12]:
        # low data rate optimization mandated for BW125 with SF11 and SF12
        DE = 1
    if sf == 6:
        # can only have implicit header with SF6
        H = 1

    Tsym = (2.0**sf)/bw
    Tpream = (Npream + 4.25)*Tsym
    payloadSymbNB = 8 + max(math.ceil((8.0*pl-4.0*sf+28+16-20*H)/(4.0*(sf-2*DE)))*(cr+4),0)
    Tpayload = payloadSymbNB * Tsym
    return Tpream + Tpayload

lamb = 0.001
T = airtime(12, 1, 20, 125) # values for experiment 4
nu = 1.0/(2*T)

for i in range(1,2001,10):
    N = i

    def f(n):
        return 2*T*N*lamb + (-2*T*(N*lamb+n*nu))*math.exp(-2*T*(N*lamb+n*nu))
    for j in range(1,10):
        print(f(j)),
    # print(" end")
    y = newton(f,2, maxiter=200)
    n = (-y/2/T - N*lamb)/nu
    print(n)
