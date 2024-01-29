use rand::{
    self,
    distributions::{Distribution, Standard},
    Rng, RngCore,
};

type Word = u64;
const NUM_WORD: usize = 8;

// little-endian
pub fn bytes_to_word(bytes: &[u8]) -> Word {
    assert!(bytes.len() < NUM_WORD);
    let mut result = 0 as Word;
    for i in 0..bytes.len() {
        if bytes[i] == 1 as u8 {
            result = result + (1 << i);
        }
    }
    result
}

pub fn word_to_bytes(word: Word) -> Vec<u8> {
    word.to_le_bytes().to_vec()
}

pub trait UniformRand: Sized {
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

impl<T> UniformRand for T
where
    Standard: Distribution<T>,
{
    #[inline]
    fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
        rng.sample(Standard)
    }
}

pub struct RngWrapper(pub rand::rngs::ThreadRng);
impl RngCore for RngWrapper {
    #[inline(always)]
    fn next_u32(&mut self) -> u32 {
        self.0.next_u32()
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }

    #[inline(always)]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.0.fill_bytes(dest)
    }

    #[inline(always)]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        self.0.try_fill_bytes(dest)
    }
}
