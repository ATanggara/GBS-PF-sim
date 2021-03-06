import numpy as np
import hafnian
from .opers import *

######## Hafnian Calculation

def sq_haf(rs):
    """
    generate cov matrix of squeezed state
    rs is array of squeezing parameters
    """
    c = np.diag(np.cosh(rs))
    s = np.diag(np.sinh(rs))
    return np.block([[c,s],
                    [s,c]])

def inter_haf(st, D):
    """
    Apply interferometer to input state
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
    compute output pattern probability from GBS with input squeezed st using hafnian
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


#### for calculation of hafnian for general state (mixed-noisy or pure)

def submtr_comp(A,n):
    """
    submatrix of complete matrix A
    """
    m = int(A.shape[0]/2)
    Du = A[:m,:m]
    Dl = A[m:,m:]
    Ou = A[:m,m:]
    Ol = A[m:,:m]
    
    Du = submtr(Du,n)
    Dl = submtr(Dl,n)
    Ou = submtr(Ou,n)
    Ol = submtr(Ol,n)
    return np.block([[Du,Ou],
                    [Ol,Dl]])

def mtr_A(s):
    """
    Construct matrix A form covariance matrix s for calculation of Hafnian function
    """
    m = int(s.shape[0]/2)
    s_Q = s + np.eye(s.shape[0])*(1/2)
    A = np.block([[np.zeros(shape=(m,m)), np.eye(m)],
                  [np.eye(m), np.zeros(shape=(m,m))]])
    s_Q_inv = np.linalg.inv(s_Q)
    A = A @ (np.eye(2*m) - s_Q_inv)
    return A


def spBS_bos(m, t, m1, m2):
    """
    Bosonic operator beamsplitter for any pair of modes
    - t transmissivity
    - m total modes of system (2 or more)
    - fm is the first mode being split
        e.g: if fm=2, then beamsplitter 
        operates on mode 2 and 3
    """
    d = np.sqrt(t) #diagonal elements
    odu = -np.sqrt(1-t) #element above the diagonal
    odl = np.sqrt(1-t) #element below the diagonal
    
    D = np.eye(m)
    D[m1-1,m1-1] = d
    D[m2-1,m2-1] = d
    D[m1-1,m2-1] = odu
    D[m2-1,m1-1] = odl
    return D


#####this only adds vacuum noise to bosonic cov matrix
def haf_therm_vac_noise(st, fm, n_bar, t, ver=False):
    """
    Add thermal noise to a mode in bosonic cov matrix
        st is covariance matrix
        fm is mode number to be mixed with termal/vacuum noise
        n_bar is the mean number of photon of thermal st
            0 for vacuum noise, -1 for no noise
        t is transmissivity of noise beamsplitter
    return: noisy cov matrix, pure cov matrix with noise mode
        if n_bar=-1, then noisy cov matrix and pure cov matrix are the same
    """
    if n_bar==-1:
        return st, st
    m = int(st.shape[0]/2)
    Du = st[:m,:m]
    Dl = st[m:,m:]
    Ou = st[:m,m:]
    Ol = st[m:,:m]
    
    ## add vacuum noise
    mean_phct = (2*n_bar +1)/2
    Du = dirsum(Du, np.eye(1)*mean_phct)
    Dl = dirsum(Dl, np.eye(1)*mean_phct)
    Ou = dirsum(Ou, np.zeros((1,1)))
    Ol = dirsum(Ol, np.zeros((1,1)))
    D = np.block([[Du,Ou],
                  [Ol,Dl]])
    if ver:
        print("D:")
        print(D)
    
    BS = spBS_bos(m+1, t, fm, m+1)
    BSl = dirsum(BS,BS)
    BSr = dirsum(BS.T,BS.T)
    D = BSl@D@BSr
    
    Du = D[:m,:m]
    Dl = D[m+1:2*m+1,m+1:2*m+1]
    Ou = D[:m,m+1:2*m+1]
    Ol = D[m+1:2*m+1,:m]
    
    return (np.block([[Du,Ou],
                      [Ol,Dl]]),
            D)

def prob_haf_gen(rs, ns, bs_arr, t, n_bar, t_noi):
    """
    compute output pattern probability using hafnian for general (pure and mixed) gaussian input state
    - rs: size (1,m) squeezing parameters
    - ns: size (1,m) output pattern
    - bs_arr: beamsplitter arrangement ([1,3,2] means beamsplit mode 1&2, 3&4, then 2&3)
    - t: transmissivity of beamsplitters
    - n_bar: amount of thermal noise.
        0 is vacuum noise, -1 no noise, otherwise thermal noise
    - t_noi: transmissivity of beamsplitters between modes and noise modes
    """
    m = rs.shape[0]
    
    #define interferometer
    D = np.eye(m)
    for j in range(bs_arr.shape[0]):
        D = spBS_bos(m=m,t=t,m1=bs_arr[j],m2=bs_arr[j]+1) @ D
        
    #calculate cov matrix
    S = sq_haf(rs)
    s = S@S.T
    s = inter_haf(st=s,D=D)
    for i in range(1,m+1): #add noise
        s,_ = haf_therm_vac_noise(st=s, fm=i, n_bar=n_bar, t=t_noi)
    
    #### calculate matrix A and A_S
    A = mtr_A(s=s)
    A_S = submtr_comp(A, ns)

    s_Q = s + np.eye(s.shape[0])*(1/2)
    det_s_Q = np.linalg.det(s_Q)
    haf = hafnian.hafnian(A_S)
    
    ##calculate constant for P(n)
    ns_fac = 1
    for i in range(ns.shape[0]):
        ns_fac = ns_fac * np.math.factorial(ns[i])
    norm_const = 1/(ns_fac*np.sqrt(det_s_Q))
    
    Pn = norm_const*np.abs(haf)
    
    return Pn


###### modify covariance matrix

def cov_mtr_reorder(s):
    """
    reorder covariance matrix form q,p,q,p to q,q,p,p
    """
    m = int(s.shape[0]/2)
    tp = type(s[0,0]) #get type of entries of s
    
    sdu = np.zeros((m,m), dtype=tp)
    sdl = np.zeros((m,m), dtype=tp)
    sou = np.ones((m,m), dtype=tp)
    sol = np.ones((m,m), dtype=tp)
    ## create upper diagonal block
    for i in range(m): #iterate over rows
        for j in range(m): #iterate over columns
            if i!=j:
                sdu[i,j] = s[2*i,2*j]
            else:
                sdu[i,j] = s[2*i,2*j]

    ## create lower diagonal block
    for i in range(m): #iterate over rows
        for j in range(m): #iterate over columns
            if i!=j:
                sdl[i,j] = s[2*i+1,2*j+1]
            else:
                sdl[i,j] = s[2*i+1,2*j+1]

    ## create upper off-diagonal block
    for i in range(m): #iterate over rows
        for j in range(m): #iterate over columns
            sou[i,j] = s[2*i,2*j+1]

    ## create lower off-diagonal block
    for i in range(m): #iterate over rows
        for j in range(m): #iterate over columns
            sol[i,j] = s[2*i+1,2*j]

    return np.block([[sdu,sou],
                  [sol,sdl]])

    
def conv_to_bos_cov(s, ver=False):
    """
    Transform quadrature convariance matrix in format (q,p,...,q,p)
        to bosonic cov matrix for GBS
    """
    m = int(s.shape[0]/2)
    T_j = (1/2)*np.array([[1, 1j],
                          [1, -1j]])
    T = T_j
    for i in range(m-1):
        T = dirsum(T_j,T)
    s_bos = T@s@T.T
    s_bos = np.real_if_close(s_bos)
    rs_bos = cov_mtr_reorder(s_bos)
    I_swap = np.block([[np.zeros((m,m)), np.eye(m)],
                       [np.eye(m), np.zeros((m,m))]])
    rs_bos = I_swap @ rs_bos
    if ver:
        print("\nbosonic cov mtr")
        print(s_bos)
        print("reordered bosonic cov mtr")
        print(rs_bos)
    return rs_bos

