"""
CV Gaussian Operations

@author: andrewtanggara
"""
import numpy as np
import hafnian

def squeezer(r):
    return np.array([[np.exp(-r), 0], [0, np.exp(r)]])

def phasor(t):
    return np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])

def dirsum(A,B):
    """
    Direct sum between matrices A and B
    """
    return np.block([[A, np.zeros((np.shape(A)[0],np.shape(B)[1]))], 
                     [np.zeros((np.shape(B)[0],np.shape(A)[1])), B]])

def msqueezer(r):
    """
    m mode squeezers
    r array or list of squeezing param for each mode
    """
    m = np.shape(r)[0]
    D = squeezer(r[0])
    for i in range(1,m):
        D = dirsum(D, squeezer(r[i]))
    return D

def twomode_squeezer(r):
    """
    Two mode squeezer
    
    """
    c = np.array([[np.cosh(r), 0],
                  [0, np.cosh(r)]])
    s = np.array([[np.sinh(r), 0],
                  [0, -np.sinh(r)]])
    return np.block([[c,s],
                     [s,c]])

def mphasor(t):
    """
    m mode phasors
    t array or list of phase shift on each mode
    """
    m = np.shape(t)[0]
    D = phasor(t[0])
    for i in range(1,m):
        D = dirsum(D, phasor(t[i]))
    return D

def beamsplitter(st, t, fm):
    """
    t in [0,1]
    fm is the first mode being split
        e.g: if fm=2, then beamsplitter 
        operates on mode 2 and 3
    st is covariance matrix
    """
    D = np.array([
        [np.sqrt(t), 0, np.sqrt(1-t), 0],
        [0, np.sqrt(t), 0, np.sqrt(1-t)],
        [-np.sqrt(1-t), 0, np.sqrt(t), 0],
        [0, -np.sqrt(1-t), 0, np.sqrt(t)]
    ])
    if (st.shape[0] > 4):
        D = dirsum(dirsum(np.eye(2*(fm-1)), D), 
                   np.eye(st.shape[0] - 2*(fm+1)))
    return D.T@st@D

def sq_haf(rs):
    c = np.diag(np.cosh(rs))
    s = np.diag(np.sinh(rs))
    return np.block([[s,c],
                    [c,s]])

def inter_haf(st, D):
    """
    D is interferometer
    st is covariance matrix
    """
    Dl = dirsum(D,D)
    Dr = dirsum(D.T,D.T)
    return (1/2)*Dl@st@Dr

def BS_bos(m, t, fm):
    """
    Bosonic operator beamsplitter
    - t transmissivity
    - m total modes of system (2 or more)
    - fm is the first mode being split
        e.g: if fm=2, then beamsplitter 
        operates on mode 2 and 3
    """
    D = np.array([[np.sqrt(t), -np.sqrt(1-t)],
                  [np.sqrt(1-t), np.sqrt(t)]])
    if (m > 2):
        D = dirsum(dirsum(np.eye((fm-1)), D), 
                   np.eye(m - (fm+1)))
    return D

def submtr(B,n):
    nidx = np.argwhere(n).flatten()
    return B[np.ix_(nidx, nidx)]

def prob_haf(rs, ns, bs_arr, t):
    """
    compute output pattern probability using hafnian
    - rs: size (1,m) squeezing parameters
    - ns: size (1,m) output pattern
    - bs_arr: beamsplitter arrangement ([1,3,2] means beamsplit mode 1&2, 3&4, then 2&3)
    - t: transmissivity of beamsplitters
    """
    m = ns.shape[0]
    #define interferometer
    D = BS_bos(m,t,1) 
    D = np.eye(m)
    for i in range(bs_arr.shape[0]):
        D = BS_bos(m, t, bs_arr[i]) @ D

    #calculate cov matrix
    s = sq_haf(rs)@sq_haf(rs).T
    s = inter_haf(s,D)

    #calculate matrix B
    B = np.diag(np.tanh(rs))
    B = D@B@D.T
    B = submtr(B, ns)

    s_Q = s + np.eye(s.shape[0])*(1/2)
    det_s_Q = np.linalg.det(s_Q)
    haf = hafnian.hafnian(B)
    Pn = (1/np.sqrt(det_s_Q)) * haf**2

    return Pn


#### EPR cov mtr
    
def cov_epr(r):
    #EPR cov matrix
    v = np.cosh(2*r)
    vI = np.array([[v, 0],
                   [0, v]])
    vZ = np.array([[np.sqrt(v**2 - 1), 0],
                   [0, -np.sqrt(v**2 - 1)]])
    s_epr = np.block([[vI, vZ],
                       [vZ, vI]])
    return s_epr

def cov_randphase_epr(r):
    #EPR cov matrix
    s_rand = cov_epr(r)
    
    #rand phase shift
    ts = np.ones(2)*np.random.rand(1)
    pi = np.ones(2)*np.pi
    s_rand = mphasor(np.divide(2*pi,ts)).T@s_rand@mphasor(np.divide(2*pi,ts))
    
    return s_rand*0.5

