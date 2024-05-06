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

########################################################## Tower Fields
class Fp(object):
    # p = 36 * x^4 + 36 * x^3 + 24 * x^2 + 6 * x + 1
    # assert(is_integer_type(p) == True)
    p = px(x)
    R = 2 ** 256
    ## Fp is constructed with 4 u64 words
    print('charactor p is {}-bits long '.format(len(bin(p)[2:])))
    R1 = R % p
    R2 = (R * R) % p
    R3 = (R1 * R2) % p
    N = R - inv_mod(p, R)
    
    alpha = Fp(-1)
    
    def ZERO():
        return Fp(0, False)
    
    def ONE():
        return Fp(1, True)
    
    def _redc(self, T):
        #assert T < (self.R*p-1)
        m = ((T & (self.R-1)) * self.N) & (self.R-1)
        t = (T + m * self.p) >> 256
        if t >= self.p:
            t -= self.p
        return t

    def __init__(self, x, redc_needed = True):
        #assert v >= 0
        if redc_needed:
            #assert v < self.R
            self.v = self._redc(x * self.R2)
            #assert self.value() == v % p
        else:
            #assert v < p
            self.v = x

    def __eq__(self, other):
        return self.v == other.v

    def __ne__(self, other):
        return self.v != other.v

    def __str__(self):
        return "%d" % (self.value())

    def __add__(self, other):
        x = self.v + other.v
        if (x >> 255) > 0 and x >= self.p:
            x -= self.p
        #assert self._redc(x) == (self.value() + other.value()) % p
        return Fp(x, False)

    def __sub__(self,other):
        x = self.v - other.v
        if x < 0:
            x += self.p
        #assert self._redc(x) == (self.value() - other.value()) % p
        return Fp(x, False)
    
    def __neg__(self):
        return Fp(self.p - self.v, False)

    def __mul__(self,other):
        return Fp(self._redc(self.v * other.v), False)

    def square(self):
        return self * self

    def double(self):
        return self + self

    def triple(self):
        return self + self + self

    def is_one(self):
        return self.v == self.R1

    def is_zero(self):
        return self.v == 0

    def value(self):
        return self._redc(self.v)

    def inverse(self):
        # Fermat
        x = Fp(self._redc(self.R3 * inv_mod_p(self.v, self.p)), False)
        #assert (self * x).value() == 1
        return x

    def additive_inverse(self):
        x = Fp(self.p, True) - self
        #assert (self + x).value() == 0
        return x

    def to_bytes(self):
        p_bytes = (self.p.bit_length() // 8) + (1 if (self.p.bit_length() % 8) > 0 else 0)
        return self.value().to_bytes(p_bytes, 'big')
    
class Fp2(object):
    ## constant beta, cubic non-residue
    beta = Fp2(Fp(1), Fp(3))
    ## Fp2 = Fp[u] / X^2 - alpha, where alpha = -1
    non_residue = Fp.alpha
    
    def ZERO():
        return Fp2(Fp.ZERO(), Fp.ZERO())
    
    def ONE():
        return Fp2(Fp.ZERO(), Fp.ONE())
    
    def __init__(self, x, y):
        """
        self.Represented as u*x + y
        """
        if type(x) == Fp:
            self.x = x
            self.y = y
        else:
            # Assumed to be integers
            self.x = Fp(x)
            self.y = Fp(y)

    def __repr__(self):
        return "(%d,%d)" % (self.x.value(), self.y.value())

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y

    def __ne__(self,other):
        return self.x != other.x or self.y != other.y

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self):
        """
        For gamma = A + iB \in gfp2
        gamma^p = A - iB
        """
        return Fp2(self.x.additive_inverse(), self.y)

    def negative_of(self):
        return Fp2(self.x.additive_inverse(), self.y.additive_inverse())

    def add(self, other):
        return Fp2((self.x + other.x), (self.y + other.y))

    def sub(self, other):
        return Fp2((self.x - other.x), (self.y - other.y))

    def double(self):
        return Fp2((self.x.double()), (self.y.double()))

    def mul(self, b):
        # assert type(a) == Fp2 and type(b) == Fp2
        # Karatsuba
        vy = (self.y * b.y)
        vx = (self.x * b.x)
        c0 = (vy - vx)
        c1 = ((self.x + self.y)*(b.x + b.y) - vy - vx)

        return Fp2(c1,c0)

    def __mul__(a,b):
        return a.mul(b)

    def __sub__(a,b):
        return a.sub(b)

    def __add__(a,b):
        return a.add(b)

    def mul_scalar(self, k):
        return Fp2((self.x * k), (self.y * k))

    # Multiply by i+3
    def mul_beta(a):
        # (beta + y)(3 + i) = 3beta + 3y - x + yi = (3x + y)i + (3y - x)
        tx = (a.x.triple()) + a.y
        ty = (a.y.triple()) - a.x
        return Fp2(tx, ty)

    def square(self):
        assert type(self.x) == Fp
        assert type(self.y) == Fp
        # Complex squaring
        t1 = self.y - self.x
        t2 = self.y + self.x
        ty = (t1 * t2)
        #ty = a.y*a.y - a.x*a.x
        tx = (self.x * self.y)
        tx = tx.double()
        return Fp2(tx, ty)

    def inverse(self):
        # Algorithm 8 from http://eprint.iacr.org/2010/354.pdf
        t = self.x.square() + self.y.square()

        inv = t.inverse()

        c_x = (self.x.additive_inverse() * inv)
        c_y = (self.y * inv)

        return Fp2(c_x, c_y)

    def exp(self, k):
        # assert is_integer_type(k)
        # assert(type(x) == Fp2)
        
        R = [Fp2(Fp(0), Fp(1)), self]
        for kb in bits_of(k):
            R[kb^1] = R[kb].mul(R[kb^1])
            assert type(R[kb]) == Fp2
            R[kb] = R[kb].square()
        return R[0]
    
# cubic extension of Fp2
class Fp6(object):
    ## gamma, which is non-quadratic residue
    gamma = Fp6(Fp2.ZERO(), Fp2.ONE(), Fp2.ZERO())
    ## Fp6 = Fp2[v] / X^3 - beta, where beta = u + 3
    non_residue = Fp2.beta
    
    def ZERO():
        return Fp6(Fp2.ZERO(), Fp2.ZERO(), Fp2.ZERO())
    
    def ONE():
        return Fp6(Fp2.ZERO(), Fp2.ZERO(), Fp2.ONE())
    
    def __init__(self, x, y, z):
        assert type(x) == Fp2 and type(y) == Fp2 and type(z) == Fp2
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __repr__(self):
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero() and self.z.is_zero()

    def is_one():
        return self.x.is_zero() and self.y.is_zero() and self.z.is_one()

    def negative_of(self):
        return Fp6(self.x.negative_of(), self.y.negative_of(), self.z.negative_of())

    def add(self, b):
        return Fp6(self.x.add(b.x),
                     self.y.add(b.y),
                     self.z.add(b.z))

    def sub(self, b):
        return Fp6(self.x.sub(b.x),
                     self.y.sub(b.y),
                     self.z.sub(b.z))

    def double(self):
        return Fp6(self.x.double(),
                     self.y.double(),
                     self.z.double())

    def mul(self, b):
        # Algorithm 13 from http://eprint.iacr.org/2010/354.pdf
        # plus some short-circuits

        if self.x.is_zero():
            if self.y.is_zero():
                return b.mul_scalar(self.z)

            t0 = (b.z * self.z)
            t1 = (b.y * self.y)

            tz = (b.x + self.y) * (self.y)
            tz -= t1
            tz = tz.mul_beta()
            tz += t0

            ty = (b.y + b.z) * (self.y + self.z)
            ty -= t0
            ty -= t1

            tx = (b.x) * (self.z)
            tx += t1

            return Fp6(tx, ty, tz)

        if b.x.is_zero():
            if b.y.is_zero():
                return self.mul_scalar(b.z)

            t0 = (self.z * b.z)
            t1 = (self.y * b.y)

            tz = (self.x + self.y) * (b.y)
            tz -= t1
            tz = tz.mul_beta()
            tz += t0

            ty = (self.y + self.z) * (b.y + b.z)
            ty -= t0
            ty -= t1

            tx = (self.x) * (b.z)
            tx += t1

            return Fp6(tx, ty, tz)

        t0 = (self.z * b.z)
        t1 = (self.y * b.y)
        t2 = (self.x * b.x)

        tz = (self.x + self.y) * (b.x + b.y)
        tz -= t1
        tz -= t2
        tz = tz.mul_beta()
        tz += t0

        ty = (self.y + self.z) * (b.y + b.z)
        ty -= t0
        ty -= t1
        ty += t2.mul_beta()

        tx = (self.x + self.z) * (b.x + b.z)
        tx -= t0
        tx += t1
        tx -= t2

        return Fp6(tx, ty, tz)

    def __mul__(a,b):
        return a.mul(b)

    def __add__(a,b):
        return a.add(b)
    def __sub__(a,b):
        return a.sub(b)

    def mul_scalar(self, k):
        assert type(k) == Fp2

        return Fp6(self.x.mul(k),
                     self.y.mul(k),
                     self.z.mul(k))

    def mul_tau(a):
        tx = a.y
        ty = a.z
        tz = a.x.mul_beta()
        return Fp6(tx, ty, tz)

    def square(self):
        # Algorithm 16 from http://eprint.iacr.org/2010/354.pdf
        ay2 = self.y.double()
        c4 = (self.z * ay2)
        c5 = self.x.square()
        c1 = c5.mul_beta() + c4
        c2 = c4 - c5
        c3 = self.z.square()
        c4 = self.x + self.z - self.y
        c5 = (ay2 * self.x)
        c4 = c4.square()
        c0 = c5.mul_beta() + c3
        c2 = c2 + c4 + c5 - c3
        n = Fp6(c2, c1, c0)
        return n

    def inverse(self):
        # Algorithm 17
        XX = self.x.square()
        YY = self.y.square()
        ZZ = self.z.square()

        XY = (self.x * self.y)
        XZ = (self.x * self.z)
        YZ = (self.y * self.z)

        A = ZZ - XY.mul_beta()
        B = XX.mul_beta() - YZ
        # There is an error in the paper for this line
        C = YY - XZ

        F = (C * self.y).mul_beta()
        F += (A * self.z)
        F += (B * self.x).mul_beta()

        F = F.inverse()

        c_x = C * F
        c_y = B * F
        c_z = A * F
        return Fp6(c_x, c_y, c_z)

## quadratic extension of Fp6
class Fp12(object):
    ## Fp12 = Fp6[w] / X^2 - gamma, where gamma = v
    non_residue = Fp6.gamma
    
    ## coefficients for one time of frobenius map
    # print(type(int((i * (((Fp.p^2) - 1) // 6)))))
    # print(type(Fp2.beta))
    beta_pi_1 = [Fp2.exp(Fp2.beta, i * ((Fp.p - 1) // 6)) for i in range(1, 6)]
    ## coefficients for two times of frobenius map
    beta_pi_2 = [Fp2.exp(Fp2.beta, i * (((Fp.p ** 2) - 1) // 6)) for i in range(1, 6)]
    
    def ZERO():
        return Fp12(Fp6.ZERO(), Fp6.ZERO())
    
    def ONE():
        return Fp12(Fp6.ZERO(), Fp6.ONE())
    
    def __init__(self, x, y = None):
        assert type(x) == Fp6
        assert type(y) == Fp6
        self.x = x
        self.y = y

    def __eq__(self,other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return "(%s,%s)" % (self.x, self.y)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self):
        return Fp12(self.x.negative_of(), self.y)

    def negative_of(self):
        return Fp12(self.x.negative_of(), self.y.negative_of())

    def frobenius(self):
        e1_x = self.x.x.conjugate_of().mul(self.beta_pi_1[4])
        e1_y = self.x.y.conjugate_of().mul(self.beta_pi_1[2])
        e1_z = self.x.z.conjugate_of().mul(self.beta_pi_1[0])

        e2_x = self.y.x.conjugate_of().mul(self.beta_pi_1[3])
        e2_y = self.y.y.conjugate_of().mul(self.beta_pi_1[1])
        e2_z = self.y.z.conjugate_of()

        return Fp12(Fp6(e1_x,e1_y,e1_z), Fp6(e2_x,e2_y,e2_z))

    def frobenius_p2(self):
        e1_x = self.x.x.mul(self.beta_pi_2[4])
        e1_y = self.x.y.mul(self.beta_pi_2[2])
        e1_z = self.x.z.mul(self.beta_pi_2[0])

        e2_x = self.y.x.mul(self.beta_pi_2[3])
        e2_y = self.y.y.mul(self.beta_pi_2[1])
        e2_z = self.y.z

        return Fp12(Fp6(e1_x,e1_y,e1_z), Fp6(e2_x,e2_y,e2_z))

    def sub(self, b):
        return Fp12(self.x - b.x, self.y - b.y)

    def mul(self, b):
        # TODO Karatsuba (algo 20)
        AXBX = self.x * b.x
        AXBY = self.x * b.y
        AYBX = self.y * b.x
        AYBY = self.y * b.y
        return Fp12(AXBY + AYBX, AYBY + AXBX.mul_tau())

    def mul_scalar(self, k):
        assert type(k) == Fp6
        return Fp12(self.x.mul(k), self.y.mul(k))

    def exp(self, k):
        # assert is_integer_type(k)

        R = [Fp12(Fp6.ZERO(), Fp6.ONE()), self]

        for kb in bits_of(k):
            R[kb^1] = R[kb].mul(R[kb^1])
            R[kb] = R[kb].square()

        return R[0]

    def square(self):
        v0 = self.x * self.y
        t = self.x.mul_tau()
        t += self.y
        ty = self.x + self.y
        ty *= t
        ty -= v0
        t = v0.mul_tau()
        ty -= t

        c_x = v0.double()
        c_y = ty

        return Fp12(c_x, c_y)

    def inverse(self):
        e = Fp12(self.x.negative_of(), self.y)

        t1 = self.x.square()
        t2 = self.y.square()
        t1 = t1.mul_tau()
        t1 = t2 - t1
        t2 = t1.inverse()

        e = e.mul_scalar(t2)
        return e
    
assert(Fp(5).square() == Fp(25))
assert(Fp(5) * (Fp(5).inverse()) == Fp.ONE())
assert(Fp2.beta * (Fp2.beta.inverse()) == Fp2.ONE())
########################################################## 


########################################################## Elliptic Curves
curve_B = Fp(3)
twist_B = Fp2(Fp(0), curve_B).mul(Fp2.beta.inverse())

###################################### traits for Elliptic Curve
def point_on_curve(point, b):
    point.force_affine()
    yy = point.y.square()
    xxx = point.x.square() * point.x
    yy -= xxx
    yy -= b
    return yy.is_zero()

def point_add(a, b):
    if a.is_infinite():
        return b
    if b.is_infinite():
        return a

    """
    http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
      Z1Z1 = a.z^2
      Z2Z2 = b.z^2
      U1 = a.x*Z2Z2
      U2 = b.x*Z1Z1
      S1 = a.y*b.z*Z2Z2
      S2 = b.y*a.z*Z1Z1
      H = U2-U1
      I = (2*H)^2
      J = H*I
      r = 2*(S2-S1)
      V = U1*I
      X3 = r^2-J-2*V
      Y3 = r*(V-X3)-2*S1*J
      Z3 = ((a.z+b.z)^2-Z1Z1-Z2Z2)*H
        """

    z1z1 = a.z.square()
    z2z2 = b.z.square()
    u1 = (z2z2 * a.x)
    u2 = (z1z1 * b.x)
    h = u2 - u1

    s1 = (a.y * b.z * z2z2)
    s2 = (b.y * a.z * z1z1)
    r = s2 - s1

    if h.is_zero() and r.is_zero():
        return a.double()

    r = r.double()
    i = h.square()
    i = i.double().double()
    j = (h * i)

    V = (u1 * i)

    c_x = (r.square() - j - V.double())
    c_y = (r * (V - c_x) - s1*j.double())

    c_z = a.z + b.z
    c_z = c_z.square()
    c_z -= z1z1
    c_z -= z2z2
    c_z *= h

    return a.__class__(c_x, c_y, c_z)

def point_double(a):
    # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    """
   
    compute A = X1^2
        compute B = Y1^2
        compute C = B^2
        compute D = 2 ((X1 + B)^2 - A - C)
        compute E = 3 A
        compute F = E^2
        compute X3 = F - 2 D
        compute Y3 = E (D - X3) - 8 C
        compute Z3 = 2 Y1 Z1
    """
    A = a.x.square()
    B = a.y.square()
    C = B.square()

    t = a.x + B
    t = t.square()

    D = (t - A - C)
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double()

    c_x = (F - D.double())
    c_y = (E * (D - c_x) - C8)
    c_z = (a.y * a.z).double()

    return a.__class__(c_x, c_y, c_z)

def point_force_affine(point):
    if point.z.is_one():
        return point

    if point.z.is_zero():
        return point

    zinv = point.z.inverse()
    zinv2 = (zinv * zinv)
    zinv3 = (zinv2 * zinv)

    point.x = point.x * zinv2
    point.y = point.y * zinv3
    point.z = point.one_element()
    
    return point

def point_scalar_mul(pt, k):
    # assert is_integer_type(k)

    if k == 0:
        return pt.__class__(pt.one_element(), pt.one_element(), pt.zero_element())

    R = pt
    e = list(reversed(to_naf(k)))[1:]
    for kb in e:
        R = R.double()
        if kb == 1:
            R = R.add(pt)
        elif kb == -1:
            R = R.add(pt.negate())
    
    return R
#####################################################################################

class G1(object):
    def __init__(self, x, y, z = Fp(1)):
        assert type(x) in [Fp]
        assert type(y) in [Fp]
        assert type(z) in [Fp]

        self.x = x
        self.y = y
        self.z = z

    def zero_element(self):
        return Fp(0)

    def one_element(self):
        return Fp(1)

    def __repr__(self):
        self.force_affine()
        return "(%d, %d)" % (self.x.value(), self.y.value())

    def is_on_curve(self):
        return point_on_curve(self, curve_B)

    def is_infinite(self):
        return self.z.is_zero()

    def add(self, b):
        return point_add(self, b)

    def double(self):
        return point_double(self)

    def force_affine(self):
        point_force_affine(self)

    def scalar_mul(self, k):
        return point_scalar_mul(self, k)
    
    def negate(self):
        return G1(self.x, -(self.y), self.z)    
    
class G2(object):
    def __init__(self, x, y, z):
        assert type(x) == Fp2 and type(y) == Fp2 and type(z) == Fp2
        self.x = x
        self.y = y
        self.z = z

    def one_element(self):
        return Fp2.ONE()

    def zero_element(self):
        return Fp2.ZERO()

    def __repr__(self):
        self.force_affine()
        return "(%s, %s)" % (self.x, self.y)

    def is_on_curve(self):
        return point_on_curve(self, twist_B)

    def is_infinite(self):
        return self.z.is_zero()

    # Add two points on the twist
    def add(a, b):
        return point_add(a,b)

    def double(a):
        return point_double(a)

    def scalar_mul(self, k):
        return point_scalar_mul(self, k)

    def force_affine(self):
        return point_force_affine(self)

    def negate(self):
        return G2(self.x, self.y.negative_of(), self.z)
        
## generator of G1: (1, -2)
g1 = G1(
    Fp(5705654538364796659083211974660230379298377357377701286821204717382476057802),
    Fp(37152507632483591484435791911354797337571904428358231246265173998464956873438),
    Fp(1)
)
assert(g1.is_on_curve() == True)

## generator of G2: 
g2 = G2(
    Fp2(
        Fp(62753311236339391648764148143567231210514074623491962919106780883433880483224),
        Fp(40543586403660527254439344329784007160991630777619071944142909139758667202651)
    ),
    Fp2(
        Fp(26582211052397732893044181121390539783060781418538774925579041820575147907402),
        Fp(12563798145120946122771505835773576822967049713945453865865029634860688995606)
    ),
    Fp2.ONE()
)
assert(g2.is_on_curve() == True)
assert(g1.scalar_mul(rx(x)).is_infinite() == True)
assert(g2.scalar_mul(rx(x)).is_infinite() == True)
#####################################################################

##################################################################### Optimal Ate Pairing
def line_func_add(r, p, q, r2):
    assert type(r) == G2
    assert type(p) == G2
    assert type(q) == G1
    assert type(r2) == Fp2

    r_t = r.z.square()
    B = p.x * r_t
    D = p.y + r.z
    D = D.square()
    D -= r2
    D -= r_t
    D *= r_t

    H = B - r.x
    I = H.square()

    E = I.double().double()

    J = H * E
    L1 = D - r.y
    L1 -= r.y

    V = r.x * E

    r_x = L1.square()
    r_x -= J
    r_x -= V.double()

    r_z = r.z + H
    r_z = r_z.square()
    r_z -= r_t
    r_z -= I

    t = V - r_x
    t *= L1
    t2 = r.y * J
    t2 = t2.double()
    r_y = t - t2

    r_out = G2(r_x, r_y, r_z)

    t = p.y + r_z
    t = t.square()
    t = t - r2
    t = t - (r_z.square())

    t2 = L1 * p.x
    t2 = t2.double()
    a = t2 - t

    c = r_z.mul_scalar(q.y).double()

    b = L1.negative_of()
    b = b.mul_scalar(q.x).double()

    return (a, b, c, r_out)

def line_func_double(r, q):
    assert type(r) == G2
    assert type(q) == G1

    # cache this?
    r_t = r.z.square()

    A = r.x.square()
    B = r.y.square()
    C = B.square()

    D = r.x + B
    D = D.square()
    D -= A
    D -= C
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double() # C*8

    r_x = F - D.double()
    r_y = E * (D - r_x) - C8

    # (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    r_z = (r.y + r.z).square() - B - r_t

    assert r_z == r.y*r.z.double()

    r_out = G2(r_x, r_y,r_z)
    #assert r_out.is_on_curve()

    a = r.x + E
    a = a.square()
    a -= (A + F + B.double().double())

    t = E * r_t
    t = t.double()
    b = t.negative_of()
    b = b.mul_scalar(q.x)

    c = r_z * r_t
    c = c.double().mul_scalar(q.y)

    return (a,b,c,r_out)

def mul_line(r, a, b, c):
    assert type(r) == Fp12
    assert type(a) == Fp2
    assert type(b) == Fp2
    assert type(c) == Fp2

    # See function fp12e_mul_line in dclxvi

    t1 = Fp6(Fp2.ZERO(), a, b)
    t2 = Fp6(Fp2.ZERO(), a, b + c)

    t1 = t1 * r.x
    t3 = r.y.mul_scalar(c)
    r.x += r.y
    r.y = t3
    r.x *= t2
    r.x -= t1
    r.x -= r.y
    r.y += t1.mul_tau()

def miller(q, p):

    import copy

    assert type(q) == G2
    assert type(p) == G1

    Q = copy.deepcopy(q)
    Q.force_affine()

    P = copy.deepcopy(p)
    P.force_affine()

    mQ = copy.deepcopy(Q)
    mQ.negate()

    f = Fp12(Fp6.ZERO(), Fp6.ONE())
    T = Q

    Qp = Q.y.square()
    
    # 6x + 2 in NAF
    naf_6xp2 = list(reversed(to_naf(6 * x + 2)))[1:]

    for naf_i in naf_6xp2:
        # Skip on first iteration?
        f = f.square()

        a, b, c, T = line_func_double(T, P)
        mul_line(f, a, b, c)

        if naf_i == 1:
            a, b, c, T = line_func_add(T, Q, P, Qp)
            mul_line(f, a, b, c)
        elif naf_i == -1:
            a, b, c, T = line_func_add(T, mQ, P, Qp)
            mul_line(f, a, b, c)

    # Q1 = pi(Q)
    Q1 = G2(
        Q.x.conjugate_of().mul(Fp12.beta_pi_1[1]),
        Q.y.conjugate_of().mul(Fp12.beta_pi_1[2]),
        Fp2.ONE())

    # Q2 = pi2(Q)
    Q2 = G2(
        Q.x.mul_scalar(Fp12.beta_pi_2[1].y),
        Q.y,
        Fp2.ONE())

    Qp = Q1.y.square()
    a, b, c, T = line_func_add(T, Q1, P, Qp)
    mul_line(f, a, b, c)

    Qp = Q2.y.square()
    a, b, c, T = line_func_add(T, Q2, P, Qp)
    mul_line(f, a, b, c)

    return f

def final_exp(inp):
    assert type(inp) == Fp12

    # Algorithm 31 from https://eprint.iacr.org/2010/354.pdf

    t1 = inp.conjugate_of()
    inv = inp.inverse()

    t1 = t1.mul(inv)
    # Now t1 = inp^(p**6-1)

    t2 = t1.frobenius_p2()
    t1 = t1.mul(t2)

    fp1 = t1.frobenius()
    fp2 = t1.frobenius_p2()
    fp3 = fp2.frobenius()

    fu1 = t1.exp(x)
    fu2 = fu1.exp(x)
    fu3 = fu2.exp(x)

    y3 = fu1.frobenius()
    fu2p = fu2.frobenius()
    fu3p = fu3.frobenius()
    y2 = fu2.frobenius_p2()

    y0 = fp1.mul(fp2)
    y0 = y0.mul(fp3)

    y1 = t1.conjugate_of()
    y5 = fu2.conjugate_of()
    y3 = y3.conjugate_of()
    y4 = fu1.mul(fu2p)
    y4 = y4.conjugate_of()

    y6 = fu3.mul(fu3p)
    y6 = y6.conjugate_of()

    t0 = y6.square()
    t0 = t0.mul(y4)
    t0 = t0.mul(y5)

    t1 = y3.mul(y5)
    t1 = t1.mul(t0)
    t0 = t0.mul(y2)
    t1 = t1.square()
    t1 = t1.mul(t0)
    t1 = t1.square()
    t0 = t1.mul(y1)
    t1 = t1.mul(y0)
    t0 = t0.square()
    t0 = t0.mul(t1)

    return t0

def optimal_ate(a, b):
    assert type(a) == G2
    assert type(b) == G1

    e = miller(a, b)
    mu = final_exp(e)

    if a.is_infinite() or b.is_infinite():
        return Fp12(Fp6.ZERO(), Fp6.ONE())

    return mu

Q, P = g2.scalar_mul(3), g1.scalar_mul(4)
mu_r = optimal_ate(Q, P)
assert(mu_r.exp(rx(x)) == Fp12.ONE())

