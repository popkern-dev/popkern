import numpy as np
from verify_stam import char_net, build_A_Abar, picks
from verify_final import vhat_fixed

def masked_const(M, P):
    """M ~ const*P entrywise ignoring 0/0. Returns (const, maxdev)."""
    mask = np.abs(P) > 1e-12
    r = M[mask]/P[mask]
    return r[0], np.abs(r - r[0]).max(), np.abs(M[~mask]).max() if (~mask).any() else 0.0

for N in [3,4,5,6,7,8,12]:
    lam, C = char_net(N)
    A, Abar = build_A_Abar(N)
    Cbar = Abar@C
    P = picks(N)
    B = {k: np.array([Cbar[i-1] for i in P[k]]).reshape(4,4) for k in (1,2,3)}
    lam2, vh = vhat_fixed(N)
    w = np.exp(2j*np.pi/N)
    g = lambda j,k: (w**(j%N))*vh[k-1]
    # PR08 grids in (a,b): rows a=0..3, cols b (patch-specific)
    # patch B (corner, [1,2]^2): a=0..3, b=0..3
    prB = np.array([[g(0,1),g(1,2),g(1,5),g(1,8)],
                    [g(0,2),g(0,3),g(0,4),g(0,13)],
                    [g(0,5),g(0,6),g(0,7),g(0,12)],
                    [g(0,8),g(0,9),g(0,10),g(0,11)]])
    # patch A (spoke-side, [1,2]x[0,1]): a=0..3, b=-1..2 ; (a,-1)=q_{j-1,{2,3,4,13}}
    prA = np.array([[g(-1,2),g(0,1),g(1,2),g(1,5)],
                    [g(-1,3),g(0,2),g(0,3),g(0,4)],
                    [g(-1,4),g(0,5),g(0,6),g(0,7)],
                    [g(-1,13),g(0,8),g(0,9),g(0,10)]])
    # patch C ([0,1]x[1,2]): a=-1..2, b=0..3 ; (-1,b)=(b,1)_{j+1}={q_{j+2,2}, q_{j+1,3}, q_{j+1,6}, q_{j+1,9}}
    prC = np.array([[g(2,2),g(1,3),g(1,6),g(1,9)],
                    [g(0,1),g(1,2),g(1,5),g(1,8)],
                    [g(0,2),g(0,3),g(0,4),g(0,13)],
                    [g(0,5),g(0,6),g(0,7),g(0,12)]])
    # Stam B_k[r][c] = refined pt (u=c, v=r) with rows r = v = Fig-1 bottom row first.
    # Mapping: PR08 (a,b) = (v_stam, u_stam). Omega2 rows v=0..3 -> a=0..3 => B2 direct vs prB.
    # Omega1: v=-? Omega1 = u in [1/2,1], v in [0,1/2]: its 4x4 rows are v=-1..2 -> a=-1..2 = patch C.
    # Omega3: u in [0,1/2]: cols u=-1..2 -> b=-1..2, rows v=0..3 -> a=0..3 = patch A.
    out = []
    for k, pr, name in ((2, prB, "O2~B"), (1, prC, "O1~C"), (3, prA, "O3~A")):
        cst, dev, z = masked_const(B[k], pr)
        out.append(f"{name}: dev={dev:.2e} zero-ok={z:.1e}")
    c2,_,_ = masked_const(B[2], prB); c1,_,_ = masked_const(B[1], prC); c3,_,_ = masked_const(B[3], prA)
    same = max(abs(c1-c2), abs(c3-c2))
    print(f"N={N}: " + "  ".join(out) + f"  |const agree|={same:.2e}")
