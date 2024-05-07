def is_py3():
    return (sys.version_info[0] == 3)

def is_integer_type(x):
    if is_py3():
        return type(x) in [int]
    else:
        return type(x) in [int,long]
    
def to_naf(x):
    z = []
    while x > 0:
        if x % 2 == 0:
            z.append(0)
        else:
            zi = 2 - (x % 4)
            x -= zi
            z.append(zi)
        x = x // 2
    return z

def bits_of(k):
    return k.bits()
    #return [int(c) for c in "{0:b}".format(k)]
    
#### traits for big integer
def gcd(x, N):
    assert(x < N)
    
    (A, B) = (N, x)
    (Ua, Ub) = (0, 1)
    (Va, Vb) = (1, 0)
    i = 0
    while B:
        q = A // B
        (A, B) = (B, A % B)
        (Ua, Ub) = (Ub, Ua + q * Ub)
        (Va, Vb) = (Vb, Va + q * Vb)
        i += 1
    r, u, v = A, Ua, Va
    
    return r, u, i % 2 == 0

def inv_mod(a, p):
    r, inv, sym = gcd(a, p)
    assert(r == 1)
    ## v * N - u * x = 1
    if sym:
        inv = p - inv
    return inv
    
def inv_mod_p(a, p):
    # Fermat
    #return (a ** (p - 2)) % p
    return pow(a, p - 2, p)

########################################################### parameter polynomials
x = ZZ['x'].gen()
## prime field, p
px = 36 * x^4 + 36 * x^3 + 24 * x^2 + 6 * x + 1
## largest prime facotr, r
rx = 36 * x^4 + 36 * x^3 + 18 * x^2 + 6 * x + 1
## trace of frobenius, t
tx = 6 * x^2 + 1
## cyclotomic group of 12-extension degree, phi_12
phi_12 = px^4 - px^2 + 1
## power final exponentiation, h
hx = (px^12 - 1) // rx
## optimal lambda in miller loop, lambda
lambdax = 6 * x + 2 + px - px^2 + px^3
## multiples of r, m
mx = lambdax // rx
#########################################################

## constant parameters
# x = 4965661367192848881
x = 6518589491078791937
# p, r, h, lamb, m = px(x), rx(x), hx(x), lambdax(x), mx(x)