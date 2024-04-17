from functools import reduce
from operator import mul
from operator import add

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

def Fp12_t_into_Fp12(e_fp12_t, F, gen):
    coef = list(e_fp12_t)
    result = []
    for i in range(len(coef)):
        result.append([(F(c) * (((gen ** 6) - F(1)) ** j) * (gen ** i)) for j, c in enumerate(coef[i].polynomial().list())])
    
    return reduce(add, sum(result, []))         

def untwist(x, y, t_x, t_y):
    return x / t_x, y / t_y

def twist(x, y, t_x, t_y):
    return x * t_x, y * t_y

def into_E12(point, beta, F, gen, t_x, t_y, E):
    x, y = list(point)[0], list(point)[1]
    x, y = into_Fp12(x, beta, F, gen), into_Fp12(y, beta, F, gen)
    return E(untwist(x, y, t_x, t_y))

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

print('\n ##################################### Curve G1: \n cofactor = {}, \n generator = {}, \n order = {} \n'.format(cofactor_E1, g_E1, r_E))

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

## twist curve E' over Fp12
Fp12_t = Fp2.extension(pol12, 'x')
Efp12_t = Efp2_t.change_ring(Fp12_t)
print('\n Twist curve E defined over Fp12: {}\n'.format(Efp12_t))
    
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
Efp12 = Efp.change_ring(Fp12)
g_E12 = into_E12(g_E2, beta, Fp, w, beta_t_x, beta_t_y, Efp12)

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
print('\n #### UNTWIST: Point of E2 \n {} \n is mapped into E12 \n {} \n successfully! \n'.format(g_E2, g_E12))

## make sure it can be twisted back
x, y = twist(list(g_E12)[0], list(g_E12)[1], beta_t_x, beta_t_y)
x, y = (into_Fp2(x, Fp, u), into_Fp2(y, Fp, u))
assert(Efp2_t(x, y) == g_E2)
print('\n #### TWIST: Point of E12 \n {} \n is mapped into E2 \n {} \n successfully! \n'.format(g_E12, Efp2_t(x, y)))

import time

## evaluation of double line divisor function
## arithmetics on fields, not on multiplicative group
def double_line(line_point, eval_point, E):
    ## lambda = 3x^2 / 2y
    (x_L, y_L) = (list(line_point)[0], list(line_point)[1])
    (x_E, y_E) = (list(eval_point)[0], list(eval_point)[1])
    alpha = (3 * x_L^2) / (2 * y_L)
    x_2L = alpha * alpha - 2 * x_L
    y_2L = -y_L - alpha * (x_2L - x_L)
    
    e_1 = y_E - y_L - alpha * (x_E - x_L)
    e_2 = x_E - x_2L

    return E(x_2L, y_2L), e_1, e_2

## evaluation of add line divisor function
## arithmetics on fields, not on multiplicative group
def add_line(line_left_point, line_right_point, eval_point, E):
    ## lambda = (y2 - y1) / (x2 - x1)
    (x_L, y_L) = (list(line_left_point)[0], list(line_left_point)[1])
    (x_R, y_R) = (list(line_right_point)[0], list(line_right_point)[1])
    (x_E, y_E) = (list(eval_point)[0], list(eval_point)[1])
    alpha = (y_L - y_R) / (x_L - x_R)
    x_LR = alpha * alpha - x_L - x_R
    y_LR = -y_L - alpha * (x_LR - x_L)
    
    e_1 = y_E - y_L - alpha * (x_E - x_L)
    e_2 = x_E - x_LR
    
    return E(x_LR, y_LR), e_1, e_2

def WeilPairing(P, Q, G, q):
    T = P
    f1 = 1
    f2 = 1
    e_bits = [int(i) for i in bin(q)[2:]]
    ## last bit cannot be evaluated, since the slope would be a vertical line
    for i in range(1, len(e_bits)):
        if (i == len(e_bits) - 1) and (e_bits[i] == 0):
            f1 = f1 * (list(Q)[0] - list(T)[0])
            T = 2 * T
            break
        T, e_1, e_2 = double_line(T, Q, G)
        f1, f2 = (f1 * f1 * e_1, f2 * f2 * e_2)
        if (i == len(e_bits) - 1) and (e_bits[i] == 1):
            f1 = f1 * (list(Q)[0] - list(T)[0])
            T = T + P
            break
        if e_bits[i] == 1:
            T, e_1, e_2 = add_line(T, P, Q, G)
            f1, f2 = (f1 * e_1, f2 * e_2)
    assert(T == G(0))
    
    return f1 / f2

G1, G2_t, G12, G12_t = (Efp, Efp2_t, Efp12, Efp12_t)
C1, C2 = (cofactor_E1, cofactor_E2_t)
## make sure they are in G1 and G2_t repectively
P, Q = (C1 * G1.random_element(), C2 * G2_t.random_element())
assert(q * P == G1(0))
assert(q * Q == G2_t(0))

###################################################
## There are three ways to obtain the Weil Pairing
###################################################

###################################### 1. map P/Q into G12 explicitly, make sure they are on the same curve(over same field)
## operations are mostly defined over Fp12
Qx = into_E12(Q, beta, Fp, w, beta_t_x, beta_t_y, G12)
assert(q * Qx == G12(0))
assert(trace_map(Qx, 12, p, G12) == G12(0))
## trival homomorphism mapping, since P is in base prime field
phi = Hom(Fp, Fp12)(Fp.gen().minpoly().roots(Fp12)[0][0])
Px = G12(phi(P.xy()[0]), phi(P.xy()[1]))
assert(Px.curve() is Qx.curve())
t1_s = time.perf_counter()
f_rP_Q_1 = WeilPairing(Px, Qx, G12, q)
f_rQ_P_1 = WeilPairing(Qx, Px, G12, q)
t1_e = time.perf_counter()
mu_1 = ((-1) ** q) * (f_rP_Q_1 / f_rQ_P_1)
assert(mu_1 ** q == Fp12(1))

###################################### 2. map P/Q into G12 implicitly
## operations are mostly defined over Fp12
assert(P.curve() is not Qx.curve())
t2_s = time.perf_counter()
f_rP_Q_2 = WeilPairing(P, Qx, G1, q)
f_rQ_P_2 = WeilPairing(Qx, P, G12, q)
t2_e = time.perf_counter()
mu_2 = ((-1) ** q) * (f_rP_Q_2 / f_rQ_P_2)
assert(mu_2 ** q == Fp12(1))
assert(f_rP_Q_1 == f_rP_Q_2)
assert(f_rQ_P_1 == f_rQ_P_2)
assert(mu_1 == mu_2)

###################################### 3. map P/Q into G12_t explicitly
## operations are mostly defined over Fp12_t / Fp2_t
## more time-consuming, since Frobenius Map based high efficient squaring/exponentiation algorithms 
## over Cyclotomic Subgroup not applied within Miller Loop, scalar multiplication(point double/add) 
## and basic field operations(evaluation of divisor function), including add/mul/inv, are what 
## Miller Loop actually do
x, y = twist(Px.xy()[0], Px.xy()[1], beta_t_x, beta_t_y)
Px = G12_t(x, y)
t3_s = time.perf_counter()
f_rP_Q_3 = WeilPairing(Px, Q, G12_t, q)
f_rQ_P_3 = WeilPairing(Q, Px, G2_t, q)
t3_e = time.perf_counter()
mu_3 = ((-1) ** q) * (f_rP_Q_3 / f_rQ_P_3)
assert(mu_3 ** q == Fp12_t(1))
mu_3 = Fp12_t_into_Fp12(mu_3, Fp, w)
assert(mu_3 ** q == Fp12(1))
assert(mu_3 == mu_2)

print('\n Time comparation: t1 = {:.2f}, t2 = {:.2f}, t3 = {:.2f}'.format(t1_e - t1_s, t2_e - t2_s, t3_e - t3_s))