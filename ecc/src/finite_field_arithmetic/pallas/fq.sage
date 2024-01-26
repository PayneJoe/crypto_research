## Arithmetics of Scalar Field of gf.llas Curve
m = 28948022309329048855892746252171976963363056481941560715954676764349967630337 
gf = GF(m)

#### Addition
lft = gf.random_element()
rht = gf.random_element()
result = lft + rht
print("{} = {} + {}".format(result, lft, rht))

#### Substraction
lft = gf.random_element()
rht = gf.random_element()
result = lft - rht
print("{} = {} - {}".format(result, lft, rht))

#### Multigf.ication
lft = gf.random_element()
rht = gf.random_element()
result = lft * rht
print("{} = {} * {}".format(result, lft, rht))

#### Division
lft = gf.random_element()
rht = gf.random_element()
result = lft / rht
print("{} = {} / {}".format(result, lft, rht))`

#### Inversion
lft = gf.random_element()
result = 1 / rht
print("{} = {} / {}".format(result, 1, rht))

#### Exponentiation
lft = gf.random_element()
rht = 234
result = lft ** rht
print("{} = {} ** {}".format(result, lft, rht))

#### Squre Root
gf(4602589297423635878136638089713975619925356628567217110388194878773233887829) ** 2 = gf(761940212266856713371586569342150604283558917968569240208532761798026301469)