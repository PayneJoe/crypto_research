## Arithmetics on Pallas Curve
p = 28948022309329048855892746252171976963363056481941560715954676764349967630337
r = 28948022309329048855892746252171976963363056481941647379679742748393362948097
Fp = GF(p)
Fr = GF(r)
b = Fp(5)
E = EllipticCurve([Fp(0), b])
g = E.gen(0)

E.order() == r

#### Addition
lft = E.random_point()
rht = E.random_point()
result = lft + rht
print("{0} \n = {1} \n + \n {2}".format(result, lft, rht))

#### ScalarMul
lft = E.random_point()
rht = 2345
result = lft * rht
print("{0} \n = {1} \n * \n {2}".format(result, lft, rht))

#### Negation
lft = E.random_point()
result = -lft
print("{0} \n = -{1}".format(result, lft))

#### Substraction
lft = E.random_point()
rht = E.random_point()
result = lft - rht
print("{0} \n = {1} \n - \n {2}".format(result, lft, rht))