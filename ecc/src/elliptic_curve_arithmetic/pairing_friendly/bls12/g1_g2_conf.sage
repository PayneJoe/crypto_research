from functools import reduce
from operator import mul

## map element of Fp2 into Fp12
def into_Fp12(e_fp2, beta, F, gen):
    a = beta.polynomial().list()
    if len(a) == 1 :
        a = a + [0]
    e = e_fp2.polynomial().list()
    if len(e) == 1:
        e = e + [0]
    return F(e[0]) + F(e[1]) * (gen ** 6 - F(a[0])) / F(a[1])

def into_Fp2(e_fp12, F, gen):
    coef = e_fp12.polynomial().list()
    zero_coeff = [1 for i in range(12) if ((len(coef) > i) and (i != 0) and (i != 6) and (F(coef[i]) == F(0)))]
    assert(reduce(mul, zero_coeff) == 1)
    
    return (F(coef[0]) + F(coef[6])) + gen * F(coef[6])

def untwist(x, y, t_x, t_y):
    return x / t_x, y / t_y

def twist(x, y, t_x, t_y):
    return x * t_x, y * t_y

def anti_trace_map(point, d, p, E):
    return d * point - trace_map(point, d, p, E)

def trace_map(point, d, p, E):
    result = point
    point_t = point
    for i in range(1, d):
        point_x, point_y = list(point_t)[0], list(point_t)[1]
        point_t = E(point_x ** p, point_y ** p)
        result = result + point_t
    return result

##########################################################################

p = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
q = 52435875175126190479447740508185965837690552500527637822603658699938581184513

A = 0
B = 4

## base prime field
Fp = GF(p)

## E1 over base prime field, map any point on Efp into the q-torsion subgroup
Efp = EllipticCurve(Fp, [A, B])
r_E = Efp.order()
cofactor_E1 = r_E // q
# g_E1 = Efp(0)
# while g_E1 == Efp(0):
#     a = Efp.random_element()
#     g_E1 = cofactor * a
g_E1 = Efp(
    2262513090815062280530798313005799329941626325687549893214867945091568948276660786250917700289878433394123885724147,
    3165530325623507257754644679249908411459467330345960501615736676710739703656949057125324800107717061311272030899084
)
assert(q * g_E1 == Efp(0))
## trace map on E1 is trival, stays on E1
assert(trace_map(g_E1, 12, p, Efp) == 12 * g_E1)

print('\n ##################################### Curve G1: \n cofactor = {}, \n generator = {}, \n order = {} \n'.format(cofactor_E1, g1, r_E))

########## Fp2 = Fp[X] / X^2 - alpha
## alpha = -1
d = 2
alpha = Fp(-1)
X = Fp['X'].gen()
pol2 = X ** d - alpha
assert(pol2.is_irreducible() == True)
Fp2 = GF(p ** d, 'u', modulus = pol2)
u = Fp2.gen()

## Fp12 = Fp2[X] / X^6 - beta
d = 6
beta = u + 1
XX = Fp2['XX'].gen()
pol12 = XX ** d - beta
assert(pol2.is_irreducible() == True)
beta_t = beta 
Efp2_t = EllipticCurve(Fp2, [A, B * beta_t])
## find the proper twisted curve, who has a q-torsion subgroup which is isomorphism with Efpk's one
if Efp2_t.order() % q != 0:
    beta_t = beta ** 5
    Efp2_t = EllipticCurve(Fp2, [A, B * beta_t])

## Fp12 = Fp[X] / X^12 - 2X^6 + 2
Fp12 = GF(p ** 12, 'w', modulus = X ** 12 - 2 * (X ** 6) + 2)
w = Fp12.gen()

## constant parameters of twist/untwist 
beta_t_x = w ** 2
beta_t_y = w ** 3

## make sure g_E2 is in the q-torsion subgroup on Efp2_t
r_E2_t = Efp2_t.order()
cofactor_E2_t = r_E2_t // q
# g_E2 = Efp2_t(0)
# while g_E2 == Efp2_t(0):
#     b = Efp2_t.random_element()
#     g_E2 = cofactor_E2_t * b
g_E2 = Efp2_t([
    [
        1265792444950586559339325656560420460408530841056393412024045461464508512562612331578200132635472221512040207420018,
        12405554917932443612178266677500354121343140278261928092817953758979290953103361135966895680930226449483176258412
    ],
    [
        3186142311182140170664472972219788815967440631281796388401764195993124196896119214281909067240924132200570679195848,
        1062539859838502367600126754068373748370820338894390252225574631210227991825937548921368149527995055326277175720251
    ],
])
assert(q * g_E2 == Efp2_t(0))
print('\n #################################### Curve G2: \n cofactor = {}, \n generator = {}, \n order = {} \n'.format(cofactor_E2_t, g_E2, r_E2_t))

## make sure g_E2 is in Fp12 first, uniform the field before untwist
g2_x, g2_y = into_Fp12(list(g_E2)[0], beta, Fp, w), into_Fp12(list(g_E2)[1], beta, Fp, w)
g2_x, g2_y = untwist(g2_x, g2_y, beta_t_x, beta_t_y)
Efp12 = Efp.change_ring(Fp12)
g_E12 = Efp12(g2_x, g2_y)

## For the convenience of do Frobenius Map within Fp2, namely (x^p, y^p)
## traditionaly need 3 steps:
## 1. untwist (x, y) to (x', y'), (x', y') = (x / beta_t_x, y / beta_t_y)
## 2. do Frobenius Map within Fp12, (x'^p, y'^p) = (x^p / beta_t_x^p, y^p / beta_t_y^p)
## 3. twist back to (x, y), (x, y) = (x'^p * beta_t_x, y'^p * beta_t_y) = (x^p / beta_t_x^{p - 1}, y^p / beta_t_y^{p - 1})
## 
## Someone may wonder why wouldn't we do Frobenius Map within Fp2 directly? 
## Since one time of Frobenius Map within Fp2, phi(P), may skip out of G2, though P belongs to G2, 
## so we must do it within the FULL EXTENSION Fp12.
##
## Caching beta_t_x^{-(p - 1)} or beta_t_y^{-(p - 1)} would be much preferable
## 
twist_frob_x = into_Fp2(1 / (beta_t_x ** (p - 1)), Fp, u)
twist_frob_y = into_Fp2(1 / (beta_t_y ** (p - 1)), Fp, u)
print('\n Twist parameters: cubic_root(beta_t)^-1 = {}, sqrt(beta_t)^-1 = {} \n'.format(beta_t_x, beta_t_y))
print('\n Twist parameters for Frobenius Map within Fp2: \n cubic_root(beta_t)^-(p - 1) = {}, \n sqrt(beta_t)^-(p - 1) = {} \n'.format(
    twist_frob_x, twist_frob_y
))

print('\n ==================================== DEBUG ====================================\n ')
## make sure g_E12 is in the zero-trace subgroup of q-torsion
assert(q * g_E12 == Efp12(0))
assert(trace_map(g_E12, 12, p, Efp12) == Efp12(0))
print('\n Point of E2 \n {} \n is mapped into E12 \n {} \n successfully! \n'.format(g_E2, g_E12))

## make sure it can be twisted back
x, y = twist(list(g_E12)[0], list(g_E12)[1], beta_t_x, beta_t_y)
x, y = (into_Fp2(x, Fp, u), into_Fp2(y, Fp, u))
assert(Efp2_t(x, y) == g_E2)