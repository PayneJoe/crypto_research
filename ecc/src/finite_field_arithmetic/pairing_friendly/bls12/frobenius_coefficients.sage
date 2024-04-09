p = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787

## Fp
Fp = GF(p)

########## Fp2 = Fp[X] / X^2 - alpha
## alpha = -1
d = 2
alpha = Fp(-1)
X = Fp['X'].gen()
pol2 = X ** d - alpha
assert(pol2.is_irreducible() == True)
Fp2 = GF(p ** d, 'u', modulus = pol2)
u = Fp2.gen()
# print('## Is multiplicative group of Fp2 cyclic? {}\n'.format(u.multiplicative_order() == (p ** 2)))

## coefficients of frobenius map over Fp2
fp2_frobenius_coefficients = [Fp(0)] * d
for i in range(1, d + 1):
    fp2_frobenius_coefficients[i % d] = (alpha ** i) ** ((p - 1)/d)
print('## Coefficients of Frobenius Map over Fp2 are:\n {0}\n'.format(fp2_frobenius_coefficients))
    
## testation
print('## Testation for coefficients of Frobenius Map over Fp2 ... \n')
a = Fp2.random_element()
for i in range(1, d + 1):
    coeff = a.polynomial().list()
    c0 = coeff[0] ** (p ** i) ## apply frobenius map in Fp
    c1 = coeff[1] ** (p ** i) ## apply frobenius map in Fp
    print('phi(a)^{0} pass? {1}'.format(i, a ** (p ** i) == Fp2([c0, c1 * fp2_frobenius_coefficients[i % d]])))
    
######## Fp6 = Fp2[X] / X^3 - beta
## beta = sqrt(alpha) + 1
d = 6
beta = u + 1
XX = Fp2['XX'].gen()
pol6 = XX ** 3 - beta
assert(pol6.is_irreducible() == True)
Fp6 = Fp2.extension(pol6, 'v')
v = Fp6.gen()

## coefficients of frobenius map over Fp6
fp6_frobenius_coefficients_c1 = [Fp2([0, 0])] * d
fp6_frobenius_coefficients_c2 = [Fp2([0, 0])] * d

cubic_c1 = [1, -u]
cubic_c2 = [1, -1]

for i in range(1, d + 1):
    M = beta ** (((p ** i) - 1)/3)
    assert(M ** 3 == cubic_c1[i % 2])
    fp6_frobenius_coefficients_c1[i % d] = M
    N = beta ** ((((p ** i) - 1) * 2)/3)
    assert(N ** 3 == cubic_c2[i % 2])
    fp6_frobenius_coefficients_c2[i % d] = N

print('\n## Coefficients of Frobenius Map over Fp6 are:\n C1 = {0}\n C2 = {1} \n'.format(
    fp6_frobenius_coefficients_c1, 
    fp6_frobenius_coefficients_c2
    )
)
    
    
## testation for coefficients of frobenius map over fp6
print('\n## Testation for coefficients of Frobenius Map over Fp6 ... \n')
a = Fp6.random_element()
# a = Fp6([Fp2([1, 2]), Fp2([1, 2]), Fp2([1, 2])])
for i in range(1, d + 1):
    coeff = list(a)
    c0 = coeff[0] ** (p ** i) ## apply frobenius map in Fp2
    c1 = coeff[1] ** (p ** i) ## apply frobenius map in Fp2
    c2 = coeff[2] ** (p ** i) ## apply frobenius map in Fp2
    print(
        'phi(a)^{0} pass? {1}'.format(
            i, 
            a ** (p ** i) == Fp6(
                [
                    c0, 
                    c1 * fp6_frobenius_coefficients_c1[i % d],
                    c2 * fp6_frobenius_coefficients_c2[i % d]
                ]
            )
        )
    )

print('\n-------------------------------------\n')

########## Fp12 = Fp6[X] / X^2 - gamma, where gamma = v
d = 12
gamma = v
XXX = Fp6['XXX'].gen()
pol12 = XXX ** 2 - gamma
assert(pol12.is_irreducible() == True)
Fp12 = Fp6.extension(pol12, 'w')
w =Fp12.gen()

## coefficients of frobenius map over Fp12
fp12_frobenius_coefficients_c1 = [Fp2([0, 0])] * 12
sextic_c1 = [1, -u]
for i in range(1, d + 1):
    M = beta ** (((p ** i) - 1)/6)
    assert(M ** 6 == sextic_c1[i % 2])
    fp12_frobenius_coefficients_c1[i % d] = M
    print('phi(a)^{0} done!'.format(i))
print('## Coefficients of Frobenius Map over Fp12 are:\n {0}\n'.format(fp12_frobenius_coefficients_c1))


## testation for coefficients of frobenius map over Fp12
print('\n## Testation for coefficients of Frobenius Map over Fp12 ... \n')
a = Fp12.random_element()
# a = Fp12([Fp6([Fp2([1, 2]), Fp2([1, 2]), Fp2([1, 2])]), Fp6([Fp2([1, 3]), Fp2([1, 3]), Fp2([1, 3])])])
for i in range(1, d + 1):
    coeff = list(a)
    c0 = coeff[0] ** (p ** i) ## apply frobenius map in Fp6
    c1 = coeff[1] ** (p ** i) ## apply frobenius map in Fp6
    print('phi(a)^{0} pass? {1}'.format(i, a ** (p ** i) == Fp12([c0, c1 * fp12_frobenius_coefficients_c1[i % d]])))