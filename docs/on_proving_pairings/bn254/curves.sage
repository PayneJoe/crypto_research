load('fields.sage')

########################################################## Elliptic Curves
curve_B = Fp(3)
twist_B = Fp2(Fp(0), curve_B).mul(beta.inverse())

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

################################################################# Testation
print('\n================ Test Module of Curves =====================\n')
t1 = (g2.is_on_curve() == True)
assert(t1 == True)
print('[Test] g2.is_on_curve()? {}\n'.format(t1))

t2 = g1.scalar_mul(rx(x)).is_infinite() == True
assert(t2 == True)
print('[Test] g1.scalar_mul(rx(x)).is_infinite()? {}\n'.format(t2))

t3 = g2.scalar_mul(rx(x)).is_infinite() == True
assert(t3 == True)
print('[Test] g2.scalar_mul(rx(x)).is_infinite()? {}\n'.format(t3))
print('\n============= End of Test Module of Curves =================\n')