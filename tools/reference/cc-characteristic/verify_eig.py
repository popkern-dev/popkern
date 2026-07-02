import numpy as np
from verify_pr08 import build_A, book_Ahat

def vhat_closed(n):
    lam = (np.cos(2*np.pi/n) + 5 + np.sqrt((np.cos(2*np.pi/n)+9)*(np.cos(2*np.pi/n)+1)))/16
    L, s = lam, np.sin(2*np.pi/n)
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
        140*L**2*(8*L-1),
        -8*L**2*(16*L-1)*(4*L+13),
        -140*L**2*(8*L-1),
        0,
        -(16*L-1)*(100*L**2+42*L-1),
        -4*L*(160*L**2+132*L-1),
        -20*L*(8*L**2+15*L+1),
        0,
        20*L*(8*L**2+15*L+1),
        4*L*(160*L**2+132*L-1)])
    vhat = ((4*L-1)*vre + 2j*s*L*(64*L-1)*vim)/13020
    return lam, vhat

print("=== eigen checks (Ahat_1 vhat = lam vhat) using corrected book matrix ===")
for n in [3,4,5,6,7,8,10,20]:
    lam, vh = vhat_closed(n)
    M = book_Ahat(n, 1, row1_variant='corrected')
    res = np.abs(M @ vh - lam*vh).max()
    # also verify v0 seed: [0, 1+conj(w), 16lam-2c-6] solves (A00-lam)v0=0
    w = np.exp(2j*np.pi/n)
    v0 = np.array([0, 1+np.conj(w), 16*lam-2*w.real-6])
    res0 = np.abs(M[:3,:3] @ v0 - lam*v0).max()
    # subdominance: full spectrum of A
    A = build_A(n)
    ev = np.sort(np.abs(np.linalg.eigvals(A)))[::-1]
    print(f"n={n}: lam={lam:.6f} |Av-lv|={res:.2e} v0seed={res0:.2e} spectrum: {ev[0]:.4f},{ev[1]:.4f},{ev[2]:.4f},{ev[3]:.4f}")
    # check proportionality of numeric eigvec of Ahat1 to closed form
    Mv, MV = np.linalg.eig(M)
    k = np.argmin(np.abs(Mv-lam))
    u = MV[:,k]; 
    # match phase/scale on component 3 (q3, real in closed form? not nec.)
    ratio = vh[2]/u[2] if abs(u[2])>1e-9 else vh[3]/u[3]
    print(f"        closed-form vs numeric eigvec: {np.abs(u*ratio - vh).max():.2e}")

print()
print("=== f(1,1,0) formula check (p.114, with errata v11) ===")
for n in [3,5,6,7]:
    lam, vh = vhat_closed(n)
    b = np.array([1/6,2/3,1/6])
    M3 = np.array([[vh[2],vh[3],vh[12]],[vh[5],vh[6],vh[11]],[vh[8],vh[9],vh[10]]])
    f110 = b @ M3 @ b
    closed = 2/29295*lam*(4*lam-1)*(139264*lam**4+170496*lam**3+112*lam**2-1476*lam-11)
    print(f"n={n}: f(1,1,0)={f110:.10f} (imag {f110.imag:.1e})  closed={closed:.10f}")
