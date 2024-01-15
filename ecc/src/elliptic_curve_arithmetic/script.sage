### Weierstrass curve with special MODULUS 2003, can be converted to Montgomery curve, and the inverse also hold 
### Weierstrass curve: y^2 = x^3 + 1132x + 278
### Montgomery curve: 899 y^2 = x^3 + 1421 x^2 + x
(a1, a3, a2, a4, a6) = (0, 0, 0, 1132, 278)
p = 2003
E = EllipticCurve([GF(p)(a1), a2, a3, a4, a6])
(b2, b4, b6, b8) = E.b_invariants()
r = E.order()
g = E.gen(0)