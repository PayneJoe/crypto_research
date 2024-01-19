// referenced from Algorithm 10.42 of "handbook of elliptic and hyperelliptic curve cryptography"
// GCD of two single precision positive integers
fn gcd(a: u16, b: u16) -> (u16, u16, u16, bool) {
    assert!(a < b);
    let (mut A, mut B) = (b, a);
    let (mut Ua, mut Ub) = (0, 1);
    let (mut Va, mut Vb) = (1, 0);
    let mut n_iter = 0_u16;
    while B != 0 {
        let q = A / B;
        (A, B) = (B, A - q * B);
        (Ua, Ub) = (Ub, Ua + q * Ub);
        (Va, Vb) = (Vb, Va + q * Vb);
        n_iter = n_iter + 1;
    }
    let (d, u, v) = (A, Ua, Va);
    (u, v, d, n_iter % 2 == 0)
}

/// positive integer cyclic group ZZ_n, reduced by a u8
#[derive(Copy, Clone, Debug)]
pub struct PositiveInteger(pub u8);

impl PositiveInteger {
    fn reduce(self, v: u16) -> u8 {
        (v % self.0 as u16) as u8
    }
}

/// functionalities of multiplicative cyclic subgroup
pub trait MultiplicativeGroup {
    fn generators(self) -> Vec<u8>;
    fn subgroup_order(self, b: u8) -> u8;
}

impl MultiplicativeGroup for PositiveInteger {
    // referenced from Corollary 4.7 of "Abstract Algebra Theory and Applications"
    fn generators(self) -> Vec<u8> {
        (1..self.0)
            .map(|v| {
                let (_, _, d, _) = gcd(v as u16, self.0 as u16);
                if d == 1 {
                    v as u8
                } else {
                    0 as u8
                }
            })
            .filter(|v| v > &0)
            .collect()
    }

    fn subgroup_order(self, b: u8) -> u8 {
        // naive reduction for any positive integer
        let norm_b = self.reduce(b as u16);
        let g = self.generators();
        assert!(g.len() > 1);
        let k: Vec<u8> = (1..self.0)
            .map(|v| {
                if self.reduce((v as u16) * (g[1] as u16)) == norm_b {
                    v
                } else {
                    0_u8
                }
            })
            .filter(|v| v > &0)
            .collect();
        assert!(k.len() >= 1);
        let (_, _, d, _) = gcd(k[0] as u16, self.0 as u16);
        self.0 / (d as u8)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generators() {
        let zz12 = PositiveInteger(12_u8);
        assert_eq!(zz12.generators(), vec![1_u8, 5_u8, 7_u8, 11_u8])
    }

    #[test]
    fn test_subgroup_order() {
        let zz60 = PositiveInteger(60_u8);
        assert_eq!(zz60.subgroup_order(56), 15_u8);
        assert_eq!(zz60.subgroup_order(55), 12_u8);
    }
}
