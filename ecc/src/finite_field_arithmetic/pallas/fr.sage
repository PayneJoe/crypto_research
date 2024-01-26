## Arithmetics of Scalar Field of Pallas Curve
m = 28948022309329048855892746252171976963363056481941647379679742748393362948097
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

#### Multiplication
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

#### Square Root
gf(1411005539847194286406933453517763875751600735672730633893962034603782311192) ** 2 == gf(23200564883514806523406480047239983027714487432999116974664880975079440510063)