import os

cur_path = os.getcwd()
load('{}/fields.sage'.format(cur_path))

########################################################## Elliptic Curves
curve_B = Fp(3)
twist_B = Fp2(Fp(0), curve_B).mul(beta.inverse())

###################################### traits for Elliptic Curve
def point_on_curve(point, b):
    p = point.force_affine()
    yy = p.y.square()
    xxx = p.x.square() * p.x
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
    u1 = z2z2.mul(a.x)
    u2 = z1z1.mul(b.x)
    h = u2.sub(u1)

    s1 = a.y.mul(b.z).mul(z2z2)
    s2 = b.y.mul(a.z).mul(z1z1)
    r = s2.sub(s1)

    if h.is_zero() and r.is_zero():
        return a.double()

    r = r.double()
    i = h.square()
    i = i.double().double()
    j = h.mul(i)

    V = u1.mul(i)

    c_x = r.square().sub(j).sub(V.double())
    c_y = r.mul(V.sub(c_x)).sub(s1.mul(j).double())

    c_z = a.z.add(b.z)
    c_z = c_z.square()
    c_z = c_z.sub(z1z1)
    c_z = c_z.sub(z2z2)
    c_z = c_z.mul(h)

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
    zinv2 = zinv.square()
    zinv3 = zinv2.mul(zinv)

    return point.__class__(point.x * zinv2, point.y * zinv3, point.one_element()) 

def point_scalar_mul(pt, k):
    assert(k >= 0)

    if k == 0:
        return pt.__class__(pt.one_element(), pt.one_element(), pt.zero_element())

    R = pt
    e = list(reversed(to_naf(k)))[1:]

    for i, kb in enumerate(e):
        R = R.double()
        # assert(R.is_on_curve())
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

    def __eq__(self, other):
        return (self.x.mul(self.z.inverse()) == other.x.mul(other.z.inverse())) and \
                (self.y.mul(self.z) == other.y) 

    def zero_element(self):
        return Fp(0)

    def one_element(self):
        return Fp(1)

    def __repr__(self):
        p = self.force_affine()
        return "(%d, %d)" % (p.x.value(), p.y.value())

    def is_on_curve(self):
        return point_on_curve(self, curve_B)

    def is_infinite(self):
        return self.z.is_zero()

    def add(self, b):
        return point_add(self, b)

    def double(self):
        return point_double(self)

    def force_affine(self):
        return point_force_affine(self)

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

    def __eq__(self, other):
        # return (self.x == other.x) and (self.y == other.y) and (self.z == other.z)
        return (self.x.mul(self.z.square().inverse()) == other.x.mul(other.z.square().inverse())) and \
        (self.y.mul(self.z.square().mul(self.z).inverse()) == other.y.mul(other.z.square().mul(other.z).inverse()))

    def one_element(self):
        return Fp2.ONE()

    def zero_element(self):
        return Fp2.ZERO()

    def __repr__(self):
        t = self.force_affine()
        return "(%s, %s)" % (t.x, t.y)

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
        
## generator of G1:
g1 = G1(
    Fp(19491323635986486980056165026003970884581302300479364565163758691834883767296),
    Fp(2503817206389621232991390790939417031444960302945150474681637705185779211401),
    Fp.ONE()
)
assert(g1.is_on_curve() == True)

## generator of G2: 
g2 = G2(
    Fp2(
        Fp(11403269112307582471523194844678173363615200121780745962919905543513926078845),
        Fp(10022529265301880767558967801827554994678953177337994173174782310334418209951)
    ),
    Fp2(
        Fp(7417909083002664933410862546938954664060641619680344911439335935535164894254),
        Fp(14403937293889182757621054345090826401263455856569175761852807173588543872656)
    ),
    Fp2.ONE()
)
assert(g2.is_on_curve() == True)

def test_curves():
    ################################################################# Testation
    print('\n================ Test Module of Curves =====================\n')

    t2 = g1.scalar_mul(rx(x)).is_infinite() == True
    assert(t2 == True)
    print('[Test] g1.scalar_mul(rx(x)).is_infinite()? {}\n'.format(t2))

    t3 = g2.scalar_mul(rx(x)).is_infinite() == True
    assert(t3 == True)
    print('[Test] g2.scalar_mul(rx(x)).is_infinite()? {}\n'.format(t3))
    # print(g2.scalar_mul(123432))
    print('\n============= End of Test Module of Curves =================\n')