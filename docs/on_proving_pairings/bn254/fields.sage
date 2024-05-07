import os

cur_path = os.getcwd()
load('{}/utils.sage'.format(cur_path))

################################################################### Tower Fields
class Fp(object):
    # p = 36 * x^4 + 36 * x^3 + 24 * x^2 + 6 * x + 1
    # assert(is_integer_type(p) == True)
    p = px(x)
    R = 2 ** 256
    ## Fp is constructed with 4 u64 words
    print('characteristic p is {}-bits long. \n'.format(len(bin(p)[2:])))
    R1 = R % p
    R2 = (R * R) % p
    R3 = (R1 * R2) % p
    N = R - inv_mod(p, R)
    
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
    

## non-quadratic residue 
alpha = Fp(-1)

class Fp2(object):
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
    
## non-cubic residue
beta = Fp2(Fp(1), Fp(3))

# cubic extension of Fp2
class Fp6(object):
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
    ## coefficients for one time of frobenius map
    beta_pi_1 = [Fp2.exp(beta, i * ((Fp.p - 1) // 6)) for i in range(1, 6)]
    ## coefficients for two times of frobenius map
    beta_pi_2 = [Fp2.exp(beta, i * (((Fp.p ** 2) - 1) // 6)) for i in range(1, 6)]
    
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
    
## non-quadratic residue
gamma = Fp6(Fp2.ZERO(), Fp2.ONE(), Fp2.ZERO())

################################################################# Testation
print('\n================ Test Module of Fields =====================\n')
t1 = (Fp(5).square() == Fp(25))
assert(t1 == True)
print('[Test] Fp(5).square() == Fp(25)? {}\n'.format(t1))

t2 = (Fp(5) * (Fp(5).inverse()) == Fp.ONE())
assert(t2 == True)
print('[Test] Fp(5) * Fp(5).inverse() == Fp.ONE()? {}\n'.format(t2))

t3 = beta * (beta.inverse()) == Fp2.ONE()
assert(t3 == True)
print('[Test] beta * beta.inverse() == Fp2.ONE()? {}\n'.format(t3))
print('\n=============== End of Test Module of Fields ==============\n')