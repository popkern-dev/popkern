import sympy as sp

lam, s = sp.symbols('lambda s', positive=True)
# c and w on the unit circle: c = cos, s = sin, w = c + i s
c = (16*lam**2 - 10*lam + 1)/(2*lam)   # book (6.12)
w = c + sp.I*s
wb = c - sp.I*s
p1,p2,p3,q1,q2,r = sp.Rational(1,64),sp.Rational(3,32),sp.Rational(9,16),sp.Rational(1,16),sp.Rational(3,8),sp.Rational(1,4)
d = 0  # i=1 frequency: delta_{i,0}=0
A00 = sp.Matrix([[0,0,0],
                 [q2*d, 2*q1*c+q2, q1*(1+wb)],
                 [r*d, r*(1+w), r]])
rows = [
 [q1*d, q1+q2*w, q2,  q1, q1*w, 0, 0],      # corrected q'_4 row
 [p2*d, 2*p1*c+p3, p2*(1+wb), p1*wb, p2, p1, 0],
 [q1*d, q1*w+q2, q2,  0, q1, q1, 0],
 [p1*d, p2*(1+w), p3, p2, p1*(1+w), p2, p1],
 [0, q2, q1*(1+wb), q1*wb, q2, q1, 0],
 [0, r, r, 0, r, r, 0],
 [0, q1, q2, q1, q1, q2, q1],
 [0, 0, r, r, 0, r, r],
 [0, q1*w, q2, q2, q1*w, q1, q1],
 [0, r*w, r, r, r*w, 0, 0]]
M = sp.zeros(13,13)
M[:3,:3] = A00
for k,row in enumerate(rows):
    for j,e in enumerate(row): M[3+k, j] = e

s2 = 1 - c**2  # s^2 substitution
def simp(e):
    e = sp.expand(e)
    e = e.subs(s**2, s2)
    return sp.simplify(sp.expand(e))

v0 = sp.Matrix([0, 1+wb, 16*lam - 2*c - 6])
# check (A00 - lam I) v0 = 0
chk = simp((A00 - lam*sp.eye(3)) @ v0)
print("v0 check:", [sp.simplify(x) for x in chk])

A10 = M[3:7, 0:3]; A11 = M[3:7, 3:7]
A20 = M[7:13, 0:3]; A21 = M[7:13, 3:7]
v1 = (A11 - lam*sp.eye(4)).solve(-A10 @ v0)
v1 = sp.Matrix([simp(x) for x in v1])
v2 = (A20 @ v0 + A21 @ v1)/lam
v2 = sp.Matrix([simp(x) for x in v2])
v = sp.Matrix.vstack(v0, v1, v2)

# Book closed form: v_book = ((4l-1) vre + 2 s l (64l-1) i vim)/13020
L = lam
vre = sp.Matrix([0,
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
vim = sp.Matrix([0,
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
vbook = ((4*L-1)*vre + 2*s*L*(64*L-1)*sp.I*vim)/13020

# The two must be proportional. Compute ratio on component 2 (q2): v[1] vs vbook[1]
ratio = simp(vbook[1]/v[1])
print("ratio (component q2):", sp.simplify(ratio))
diffs = []
for k in range(13):
    dk = simp(vbook[k] - ratio*v[k])
    diffs.append(sp.simplify(dk))
for k,dk in enumerate(diffs):
    print(f"k={k+1}: diff = {dk}")
