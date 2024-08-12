pub trait BitVectorize {
    fn to_le_bits(self) -> Vec<u8>;
}

impl BitVectorize for u128 {
    fn to_le_bits(self) -> Vec<u8> {
        let mut bits = vec![0_u8; 128];
        let mut v = self.clone();
        let mut most_significant = 0;
        (0..128).for_each(|i| {
            let bit = (v & 1) as u8;
            bits[i] = bit;
            if bit == 1 {
                most_significant = i;
            }
            v = v >> 1;
        });
        bits = bits[0..(most_significant + 1)].to_vec();
        bits.reverse();
        bits
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_to_le_bits() {
        let a = 17_u128;
        let hint = vec![1_u8, 0_u8, 0_u8, 0_u8, 1_u8];
        assert!(a.to_le_bits() == hint);
    }
}
