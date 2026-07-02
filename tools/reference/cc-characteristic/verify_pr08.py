import numpy as np

def build_A(n, alpha=None, beta=None, gamma=None):
    """13n x 13n CC subdivision matrix on PR08's labeled net.
    Global index of q_{j,k} = 13*j + (k-1). Grid: sector j coords (a,b),
    P: k -> (a,b): 1:(0,0) 2:(1,0) 3:(1,1) 4:(1,2) 5:(2,0) 6:(2,1) 7:(2,2)
       8:(3,0) 9:(3,1) 10:(3,2) 11:(3,3) 12:(2,3) 13:(1,3)
    (0,b) of sector j = q_{j+1, {2,5,8}[b-1]}  (b=1,2,3)."""
    if alpha is None:
        alpha = 1 - 7/(4*n); beta = 3/(2*n); gamma = 1/(4*n)
    A = np.zeros((13*n, 13*n))
    def gi(j,k): return 13*(j % n) + k - 1
    for j in range(n):
        jm, jp = (j-1)%n, (j+1)%n
        def add(k_new, terms):
            for (jj,kk),w in terms.items(): A[gi(j,k_new), gi(jj,kk)] += w
        # k=1 EV: alpha*q1 + beta/n * sum spoke + gamma/n * sum diag
        t = {(j,1): alpha}
        for jj in range(n):
            t[(jj,2)] = t.get((jj,2),0) + beta/n
            t[(jj,3)] = t.get((jj,3),0) + gamma/n
        add(1, t)
        # k=2 edge point spoke j: 3/8(q1+qj2)+1/16(qj3+q_{j+1,2}+q_{j-1,3}+q_{j-1,2})
        add(2, {(j,1):3/8,(j,2):3/8,(j,3):1/16,(jp,2):1/16,(jm,3):1/16,(jm,2):1/16})
        # k=3 face point quad j: 1/4(q1+qj2+qj3+q_{j+1,2})
        add(3, {(j,1):1/4,(j,2):1/4,(j,3):1/4,(jp,2):1/4})
        # k=4 edge point of edge (0,1)-(1,1): 3/8(q_{j+1,2}+qj3)+1/16(q1+qj2+q_{j+1,5}+qj4)
        add(4, {(jp,2):3/8,(j,3):3/8,(j,1):1/16,(j,2):1/16,(jp,5):1/16,(j,4):1/16})
        # k=5 vertex point at (1,0)=qj2: 9/16 qj2 + 3/32(q1,qj5,qj3,q_{j-1,3}) + 1/64(q_{j+1,2},q_{j-1,2},qj6,q_{j-1,4})
        add(5, {(j,2):9/16,(j,1):3/32,(j,5):3/32,(j,3):3/32,(jm,3):3/32,
                (jp,2):1/64,(jm,2):1/64,(j,6):1/64,(jm,4):1/64})
        # k=6 edge point (1,0)-(1,1): 3/8(qj2+qj3)+1/16(q1+q_{j+1,2}+qj5+qj6)
        add(6, {(j,2):3/8,(j,3):3/8,(j,1):1/16,(jp,2):1/16,(j,5):1/16,(j,6):1/16})
        # k=7 vertex point at (1,1)=qj3: 9/16 qj3 + 3/32(qj2,q_{j+1,2},qj4,qj6) + 1/64(q1,qj5,qj7,q_{j+1,5})
        add(7, {(j,3):9/16,(j,2):3/32,(jp,2):3/32,(j,4):3/32,(j,6):3/32,
                (j,1):1/64,(j,5):1/64,(j,7):1/64,(jp,5):1/64})
        # k=8 edge point (1,0)-(2,0): 3/8(qj2+qj5)+1/16(qj3+qj6+q_{j-1,3}+q_{j-1,4})
        add(8, {(j,2):3/8,(j,5):3/8,(j,3):1/16,(j,6):1/16,(jm,3):1/16,(jm,4):1/16})
        # k=9 face point {(1,0),(2,0),(2,1),(1,1)}: 1/4(qj2+qj5+qj6+qj3)
        add(9, {(j,2):1/4,(j,5):1/4,(j,6):1/4,(j,3):1/4})
        # k=10 edge point (1,1)-(2,1): 3/8(qj3+qj6)+1/16(qj2+qj5+qj4+qj7)
        add(10,{(j,3):3/8,(j,6):3/8,(j,2):1/16,(j,5):1/16,(j,4):1/16,(j,7):1/16})
        # k=11 face point {(1,1),(2,1),(2,2),(1,2)}: 1/4(qj3+qj6+qj7+qj4)
        add(11,{(j,3):1/4,(j,6):1/4,(j,7):1/4,(j,4):1/4})
        # k=12 edge point (1,1)-(1,2): 3/8(qj3+qj4)+1/16(q_{j+1,2}+q_{j+1,5}+qj6+qj7)
        add(12,{(j,3):3/8,(j,4):3/8,(jp,2):1/16,(jp,5):1/16,(j,6):1/16,(j,7):1/16})
        # k=13 face point {(0,1),(1,1),(1,2),(0,2)}: 1/4(q_{j+1,2}+qj3+qj4+q_{j+1,5})
        add(13,{(jp,2):1/4,(j,3):1/4,(j,4):1/4,(jp,5):1/4})
    return A

def book_Ahat(n, i, alpha=None, beta=None, gamma=None, row1_variant='printed'):
    """PR08 p.111-112 printed \hat A_i (13x13), i-th frequency block."""
    if alpha is None:
        alpha = 1 - 7/(4*n); beta = 3/(2*n); gamma = 1/(4*n)
    w = np.exp(2j*np.pi*i/n); wb = np.conj(w); c = w.real
    d = 1.0 if i % n == 0 else 0.0
    p1,p2,p3,q1,q2,r = 1/64,3/32,9/16,1/16,3/8,1/4
    M = np.zeros((13,13), complex)
    M[0,:3] = [alpha*d, beta*d, gamma*d]
    M[1,:3] = [q2*d, 2*q1*c+q2, q1*(1+wb)]
    M[2,:3] = [r*d, r*(1+w), r]
    rows = [
      [q2*d, q1+q2*w, q1,      q1,   q1*w,      0,  0],   # printed row for q'_4
      [p2*d, 2*p1*c+p3, p2*(1+wb), p1*wb, p2,   p1, 0],
      [q1*d, q1*w+q2, q2,      0,    q1,        q1, 0],
      [p1*d, p2*(1+w), p3,     p2,   p1*(1+w),  p2, p1],
      [0,    q2, q1*(1+wb),    q1*wb,q2,        q1, 0],
      [0,    r,  r,            0,    r,         r,  0],
      [0,    q1, q2,           q1,   q1,        q2, q1],
      [0,    0,  r,            r,    0,         r,  r],
      [0,    q1*w, q2,         q2,   q1*w,      q1, q1],
      [0,    r*w, r,           r,    r*w,       0,  0]]
    if row1_variant == 'corrected':
        rows[0] = [q1*d, q1+q2*w, q2, q1, q1*w, 0, 0]
    for k,row in enumerate(rows):
        M[3+k,:7] = row
    return M

for n in [3,5,6,7]:
    A = build_A(n)
    assert np.allclose(A.sum(axis=1), 1), "rows must sum to 1"
    # DFT: Ahat_i = sum_m A_m w^{im}, A_m = block mapping Q_{j+m} contributions
    A0j = [A[0:13, 13*m:13*m+13] for m in range(n)]  # row block j=0: A_m = A[0, block m] since Q'_0 = sum_m A[0,m] Q_m and A_m multiplies Q_{0+m}
    for i in [0,1,2]:
        Ahat = sum(A0j[m]*np.exp(2j*np.pi*i*m/n) for m in range(n))
        Mp = book_Ahat(n,i,row1_variant='printed')
        Mc = book_Ahat(n,i,row1_variant='corrected')
        dp = np.abs(Ahat[:, :7] - Mp[:, :7]).max()
        dc = np.abs(Ahat[:, :7] - Mc[:, :7]).max()
        zero_tail = np.abs(Ahat[:,7:]).max()
        print(f"n={n} i={i}: |Ahat-printed|={dp:.2e}  |Ahat-corrected|={dc:.2e}  tail(cols8-13)={zero_tail:.2e}")
