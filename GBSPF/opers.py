"""
CV Gaussian Operations

@author: andrewtanggara
"""
import numpy as np
import hafnian

def squeezer(r):
    return np.array([[np.exp(-r), 0], 
                      [0, np.exp(r)]])

def phasor(t):
    return np.array([[np.cos(t), np.sin(t)], 
                      [-np.sin(t), np.cos(t)]])

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
    """
    Submatrix of B based on n
        e.g: for n = (0,1,1,0), we have submatrix:
        [[B[1,1],B[1,2]], 
        [B[2,1],B[2,2]]
    """
    nidx = np.argwhere(n).flatten()
    return B[np.ix_(nidx, nidx)]

def prob_haf(rs, ns, bs_arr, t, n_bar, t_noi):
    """
    compute output pattern probability using hafnian
    - rs: size (1,m) squeezing parameters
    - ns: size (1,m) output pattern
    - bs_arr: beamsplitter arrangement ([1,3,2] means beamsplit mode 1&2, 3&4, then 2&3)
    - t: transmissivity of beamsplitters
    - n_bar: amount of thermal noise.
        0 is vacuum noise, -1 no noise, otherwise thermal noise
    - t_noi: transmissivity of beamsplitters between modes and noise modes
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
    for i in range(m):
        s,_ = therm_vac_noise(st=s, fm=i, n_bar=n_bar, t=t_noi)

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
    s_rand = mphasor(np.multiply(2*pi,ts)).T@s_rand@mphasor(np.multiply(2*pi,ts))
    
    return s_rand*0.5


######## Noise

def spbeamsplitter(st, t, m1, m2):
    """
    m1 is mode # of first input to bs
    m2 is mode # of second input to bs
    """
    m1 = m1-1 #match indexing that starts from 0
    m2 = m2-1
    D = np.array([[np.sqrt(t), 0],
                  [0, np.sqrt(t)]]) #splitter block in matrix diagonal
    ODu = np.array([[np.sqrt(1-t), 0],
                  [0, np.sqrt(1-t)]]) #splitter block in upper section of matrix
    ODl = np.array([[-np.sqrt(1-t), 0],
                  [0, -np.sqrt(1-t)]]) #splitter block in lower section of matrix
    BS = np.block([[D, ODu],
                   [ODl, D]]) #beamsplitter for 2 modes
    if (st.shape[0] > 4): #if more than 2 modes, build BS matrix on chosen modes m1 and m2
        BS = []
        for j in range(int(st.shape[0]/2)):
            row = []
            for k in range(int(st.shape[0]/2)):
                if ((j==m1) and (k==m1)) or ((j==m2) and (k==m2)):
                    row.append(D)
                elif (j==m1) and (k==m2):
                    row.append(ODu)
                elif (j==m2) and (k==m1):
                    row.append(ODl)
                elif (j==k):
                    row.append(np.eye(2))
                else:
                    row.append(np.zeros(shape=(2,2)))
            BS.append(row)
        BS = np.block(BS)
    return BS.T@st@BS

def therm(st, m, n_bar):
    """
    Add thermal noise to a mode in cov matrix.
        st is covariance matrix
        m is mode number of vacuum to be converted to thermal st.
        n_bar is the mean number of photon of thermal st.
    """
    thst = np.array([[(2*n_bar), 0],
                     [0, (2*n_bar)]])
    D = thst
    if (st.shape[0] > 2):
        D = dirsum(dirsum(np.zeros((2*(m-1), 2*(m-1))), thst), 
                   np.zeros((st.shape[0] - 2*(m), st.shape[0] - 2*(m))))
    return st + D

def therm_vac_noise(st, fm, n_bar, t):
    """
    Add thermal noise to a mode in cov matrix/
        st is covariance matrix
        fm is mode number to be mixed with termal/vacuum noise
        n_bar is the mean number of photon of thermal st
            0 for vacuum noise, -1 for no noise
        t is transmissivity of noise beamsplitter
    return: noisy cov matrix, pure cov matrix with noise mode
        if n_bar=-1, then noisy cov amtrix and pure cov matrix are the same
    """
    if n_bar==-1:
        return st, st
    D = dirsum(st, np.eye(2)) #add vacuum noise mode to cov matrix
    m = int(D.shape[0]/2) #get mode number of vacuum noise mode
    D = therm(D, m, n_bar) #add thermal noise the vacuum noise mode
    D = spbeamsplitter(D, t, fm, m) #mix noise with specified mode
    return D[:(m*2)-2, :(m*2)-2], D

