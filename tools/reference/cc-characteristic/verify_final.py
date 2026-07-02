import numpy as np
from verify_pr08 import build_A, book_Ahat

def vhat_fixed(n):
    cn = np.cos(2*np.pi/n); s = np.sin(2*np.pi/n)
    lam = (cn + 5 + np.sqrt((cn+9)*(cn+1)))/16
    L = lam
    vre = np.array([0,
        4*L**2*(64*L-1)*(32*L-1)*(16*L-1)*(4*L-1),
        8*L**2*(64*L-1)*(32*L-1)*(16*L-1),
        4*L**2*(64*L-1)*(928*L**2+228*L-31),
        8*L**2*(64*L-1)*(16*L-1)*(4*L-1)*(4*L+13),
        4*L**2*(64*L-1)*(928*L**2+228*L-31),
        80*L**2*(1280*L**3+2128*L**2-56*L-13),
        (64*L-1)*(16*L-1)*(4*L-1)*(100*L**2+42*L-1),
        4*L*(64*L-1)*(640*L**3+688*L**2-82*L-1),
        20*L*(2048*L**4+11040*L**3+812*L**2-165*L-1),
        40*L*(5248*L**3+1568*L**2-133*L-5),
        20*L*(2048*L**4+11040*L**3+812*L**2-165*L-1),
        4*L*(64*L-1)*(640*L**3+688*L**2-82*L-1)])
    vim = np.array([0,
        -4*L**2*(32*L-1)*(16*L-1),
        0,
        84*L**2*(8*L+1),          # ERRATUM fix (book: 140 L^2 (8L-1))
        -8*L**2*(16*L-1)*(4*L+13),
        -84*L**2*(8*L+1),         # ERRATUM fix
        0,
        -(16*L-1)*(100*L**2+42*L-1),
        -4*L*(160*L**2+132*L-1),
        -20*L*(8*L**2+15*L+1),
        0,
        20*L*(8*L**2+15*L+1),
        4*L*(160*L**2+132*L-1)])
    return lam, ((4*L-1)*vre + 2j*s*L*(64*L-1)*vim)/13020

def build_A_distributed(n):
    """Same as build_A but every stencil weight on the EV is spread as w/n over all copies q_{j',1}."""
    A = build_A(n)
    for row in range(13*n):
        wev = sum(A[row, 13*j] for j in range(n))
        for j in range(n): A[row, 13*j] = wev/n
    return A

print("=== residual of FIXED closed-form eigenvector, and true spectrum (distributed center) ===")
for n in [3,4,5,6,7,8,10,20,50]:
    lam, vh = vhat_fixed(n)
    M = book_Ahat(n,1,row1_variant='corrected')
    res = np.abs(M@vh - lam*vh).max() / np.abs(vh).max()
    Ad = build_A_distributed(n)
    assert np.allclose(Ad.sum(axis=1),1)
    ev = np.linalg.eigvals(Ad)
    ev = ev[np.argsort(-np.abs(ev))]
    lam_num = np.abs(ev[1])
    print(f"n={n:3d}: lam={lam:.12f} lam_num={lam_num:.12f} |Av-lv|/|v|={res:.2e} "
          f"top6=|{'|'.join(f'{abs(e):.4f}' for e in ev[:6])}|")

# subdominant eigenvector of full matrix vs W1^{-1} vhat (v_j = w^j vhat / n)
print()
print("=== full-matrix eigenvector consistency: v_{j,k} = w_n^j vhat_k (up to scale) ===")
for n in [3,5,7]:
    lam, vh = vhat_fixed(n)
    Ad = build_A_distributed(n)
    v = np.concatenate([np.exp(2j*np.pi*j/n)*vh for j in range(n)])
    print(f"n={n}: |A v - lam v| = {np.abs(Ad@v - lam*v).max():.2e}")

# f(1,1,0) sign/normalization with fixed vector
print()
b = np.array([1/6,2/3,1/6])
for n in [3,5,6,7]:
    lam, vh = vhat_fixed(n)
    M3 = np.array([[vh[2],vh[3],vh[12]],[vh[5],vh[6],vh[11]],[vh[8],vh[9],vh[10]]])
    f110 = b@M3@b
    closed = 2/29295*lam*(4*lam-1)*(139264*lam**4+170496*lam**3+112*lam**2-1476*lam-11)
    print(f"n={n}: f(1,1,0) = {f110:.12f}  closed formula = {closed:.12f}")
