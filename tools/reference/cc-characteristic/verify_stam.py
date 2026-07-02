import numpy as np

d, e, f = 3/8, 1/16, 1/4
a, b, c = 9/16, 3/32, 1/64

def build_S(N):
    """Stam S: (2N+1)x(2N+1). Vertex 1=EV; ring: edge 2j (j=1..N), face 2j+1."""
    K1 = 2*N+1
    S = np.zeros((K1,K1))
    aN, bN, cN = 1-7/(4*N), 3/(2*N**2), 1/(4*N**2)
    S[0,0] = aN
    for j in range(1,N+1):
        S[0, 2*j-1] = bN   # edge vertex 2j -> index 2j-1 (0-based)
        S[0, 2*j % K1 if 2*j<K1 else 0] = S[0,2*j] if False else S[0,2*j]  # placeholder
    # redo cleanly with helper
    def idx(v):  # vertex label (1..2N+1) -> 0-based, ring cyclic on 2..2N+1
        v = ((v - 2) % (2*N)) + 2
        return v-1
    S[:] = 0
    S[0,0] = aN
    for j in range(1,N+1):
        S[0, idx(2*j)] += bN
        S[0, idx(2*j+1)] += cN
    for j in range(1,N+1):
        r = idx(2*j)  # new edge point row
        S[r,0] += d; S[r, idx(2*j)] += d
        for v in (2*j-2, 2*j-1, 2*j+1, 2*j+2):
            S[r, idx(v)] += e
        r = idx(2*j+1)  # new face point row
        for v in (1,):
            S[r,0] += f
        for v in (2*j, 2*j+1, 2*j+2):
            S[r, idx(v)] += f
    return S

def stam_blocks(N):
    """S11 (7 x 2N+1), S12 (7x7), S21 (9 x 2N+1), S22 (9x7) per Appendix A."""
    K1 = 2*N+1
    S11 = np.zeros((7,K1)); S12 = np.zeros((7,7))
    S21 = np.zeros((9,K1)); S22 = np.zeros((9,7))
    # S12 rows (cols 2N+2..2N+8):
    S12[0] = [c,b,c,0,b,c,0]
    S12[1] = [0,e,e,0,0,0,0]
    S12[2] = [0,c,b,c,0,0,0]
    S12[3] = [0,0,e,e,0,0,0]
    S12[4] = [0,0,0,0,e,e,0]
    S12[5] = [0,0,0,0,c,b,c]
    S12[6] = [0,0,0,0,0,e,e]
    # S11 rows, first 8 cols (rest zero); N=3 special-cased below
    S11rows = [
      [c,0,0,b,a,b,0,0],
      [e,0,0,e,d,d,0,0],
      [b,0,0,c,b,a,b,c],
      [e,0,0,0,0,d,d,e],
      [e,0,0,d,d,e,0,0],
      [b,c,b,a,b,c,0,0],
      [e,e,d,d,0,0,0,0]]
    for i,row in enumerate(S11rows):
        for j,x in enumerate(row):
            if N==3 and j==7:  # col 8 folds onto col 2 (c8==c2)
                S11[i,1] += x
            else:
                S11[i,j] += x
    if N==3:
        # paper: second column of S11 equals (0,0,c,e,0,c,e)^T -- check against fold
        pass
    S21rows = [
      [0,0,0,0,f,0,0,0],
      [0,0,0,0,d,e,0,0],
      [0,0,0,0,f,f,0,0],
      [0,0,0,0,e,d,e,0],
      [0,0,0,0,0,f,f,0],
      [0,0,0,e,d,0,0,0],
      [0,0,0,f,f,0,0,0],
      [0,0,e,d,e,0,0,0],
      [0,0,f,f,0,0,0,0]]
    for i,row in enumerate(S21rows):
        for j,x in enumerate(row):
            if N==3 and j==7: S21[i,1] += x
            else: S21[i,j] += x
    S22[0] = [f,f,0,0,f,0,0]
    S22[1] = [e,d,e,0,e,0,0]
    S22[2] = [0,f,f,0,0,0,0]
    S22[3] = [0,e,d,e,0,0,0]
    S22[4] = [0,0,f,f,0,0,0]
    S22[5] = [e,e,0,0,d,e,0]
    S22[6] = [0,0,0,0,f,f,0]
    S22[7] = [0,0,0,0,e,d,e]
    S22[8] = [0,0,0,0,0,f,f]
    return S11,S12,S21,S22

def build_A_Abar(N):
    K1 = 2*N+1; K = 2*N+8; M = 2*N+17
    S = build_S(N); S11,S12,S21,S22 = stam_blocks(N)
    A = np.zeros((K,K)); Abar = np.zeros((M,K))
    A[:K1,:K1] = S
    A[K1:,:K1] = S11; A[K1:,K1:] = S12
    Abar[:K,:] = A
    Abar[K:,:K1] = S21; Abar[K:,K1:] = S22
    return A, Abar

def picks(N):
    P = {}
    P[1] = [8,7,2*N+5,2*N+13, 1,6,2*N+4,2*N+12, 4,5,2*N+3,2*N+11, 2*N+7,2*N+6,2*N+2,2*N+10]
    P[2] = [1,6,2*N+4,2*N+12, 4,5,2*N+3,2*N+11, 2*N+7,2*N+6,2*N+2,2*N+10, 2*N+16,2*N+15,2*N+14,2*N+9]
    P[3] = [2,1,6,2*N+4, 3,4,5,2*N+3, 2*N+8,2*N+7,2*N+6,2*N+2, 2*N+17,2*N+16,2*N+15,2*N+14]
    if N==3: P[1][0] = 2  # q^1_1 = 2 instead of 8 for N=3
    return P

def char_net(N):
    """K-vector complex: frequency-1 subdominant eigenvector (u4 + i u5) extended to K=2N+8.
       Outer 7 get filled by requiring A-eigenvector: solve (S12 - lam I) x = -S11 v_ring."""
    cN = np.cos(2*np.pi/N)
    lam = (cN+5+np.sqrt((cN+9)*(cN+1)))/16
    K1 = 2*N+1
    v = np.zeros(K1, complex)
    w = np.exp(2j*np.pi/N)
    for j in range(1,N+1):
        v[2*j-1]   = (4*lam-1)*w**(j-1)   # edge vertex 2j
        v[2*j % K1 if 2*j < K1 else 0] = v[2*j] if False else v[2*j]
    # redo cleanly
    v[:] = 0
    for j in range(1,N+1):
        v[(2*j)-1] = (4*lam-1)*w**(j-1)      # 0-based index of vertex 2j
        v[(2*j+1)-1 if 2*j+1<=K1 else 1] = w**(j-1)+w**j   # vertex 2j+1; for j=N: 2N+1 fine
    S = build_S(N)
    assert np.abs(S@v - lam*v).max() < 1e-12, np.abs(S@v-lam*v).max()
    S11,S12,S21,S22 = stam_blocks(N)
    x = np.linalg.solve(S12 - lam*np.eye(7), -S11@v)
    return lam, np.concatenate([v,x])

print("=== Stam S row sums, eigenvalue formulas ===")
for N in [3,4,5,6,7,8]:
    S = build_S(N)
    A, Abar = build_A_Abar(N)
    assert np.allclose(S.sum(1),1) and np.allclose(A.sum(1),1) and np.allclose(Abar.sum(1),1)
    ev = np.sort(np.abs(np.linalg.eigvals(S)))[::-1]
    # formulas
    mu23 = [(-7+3*N + s*np.sqrt(49-30*N+5*N*N+0j))/(8*N) for s in (-1,1)]
    lams = []
    for l in range(1,N):
        cl = np.cos(2*np.pi*l/N)
        for s in (-1,1):
            lams.append((5+cl + s*np.cos(np.pi*l/N)*np.sqrt(18+2*cl))/16)
    pred = np.sort(np.abs(np.array([1]+[m.real for m in mu23]+lams)))[::-1]
    print(f"N={N}: |spec(S)-formula|={np.abs(ev-pred).max():.2e}", end="  ")
    lam, C = char_net(N)
    print(f"lam={lam:.9f} (subdominant: {ev[1]:.9f})")

print()
print("=== characteristic-map invariance: chi(A C) = lam * chi(C), and picking consistency ===")
# Also: after Abar, pick the three 4x4 grids and compare Omega2 corner-grid with PR08 patch B.
from verify_final import vhat_fixed
for N in [3,5,6,7]:
    lam, C = char_net(N)   # K-vector complex char net
    A, Abar = build_A_Abar(N)
    # eigen check on extended net
    print(f"N={N}: |A C - lam C| = {np.abs(A@C - lam*C).max():.2e}", end="  ")
    Cbar = Abar@C
    P = picks(N)
    B = {k: np.array([Cbar[i-1] for i in P[k]]).reshape(4,4) for k in (1,2,3)}
    # C2 join check: Omega1 grid rows 1..3 (v-rows) == Omega2 grid rows 0..2 (shifted by one v-row)
    d12 = np.abs(B[1][1:,:] - B[2][:3,:]).max()
    d23 = np.abs(B[2][:, :3] - B[3][:, 1:]).max()
    print(f"grid-overlap O1/O2: {d12:.2e}  O2/O3: {d23:.2e}", end="  ")
    # eigen-scaling of picked grids: pick at level 2 = lam * pick at level 1
    Cbar2 = Abar@(A@C)
    B2 = {k: np.array([Cbar2[i-1] for i in P[k]]).reshape(4,4) for k in (1,2,3)}
    sc = max(np.abs(B2[k] - lam*B[k]).max() for k in (1,2,3))
    print(f"|pick(A C) - lam pick(C)| = {sc:.2e}")

print()
print("=== cross-check with PR08 net: Stam corner grid Omega2 vs PR08 patch B grid ===")
for N in [3,5,6,7]:
    lam, C = char_net(N)
    A, Abar = build_A_Abar(N)
    Cbar = Abar@C
    P = picks(N)
    B2 = np.array([Cbar[i-1] for i in P[2]]).reshape(4,4)
    # PR08: vhat (13/sector, sector j = w^j vhat). Patch B grid (a=0..3 rows, b=0..3 cols):
    lam2, vh = vhat_fixed(N)
    w = np.exp(2j*np.pi/N)
    g = lambda j,k: (w**j)*vh[k-1]
    prB = np.array([
      [g(0,1),  g(1,2), g(1,5), g(1,8)],
      [g(0,2),  g(0,3), g(0,4), g(0,13)],
      [g(0,5),  g(0,6), g(0,7), g(0,12)],
      [g(0,8),  g(0,9), g(0,10),g(0,11)]])
    # Stam B2 is the corner patch of the L at level 1 (scale lam^?); PR08 patch B at ring 0.
    # They should be related by a single complex scale (rotation+scaling), possibly transpose.
    for name, M in [("as-is", B2), ("transpose", B2.T), ("conj", np.conj(B2)), ("conjT", np.conj(B2).T)]:
        ratios = M/prB
        if np.allclose(ratios, ratios[0,0], rtol=1e-8, atol=1e-10):
            print(f"N={N}: Omega2 {name} = (complex const {ratios[0,0]:.6f}) * PR08-patchB  MATCH")
            break
    else:
        print(f"N={N}: no simple match; ratio matrix (as-is):")
        print(np.round(B2/prB, 6))
