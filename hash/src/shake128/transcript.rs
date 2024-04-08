// for the purpose of research we use Shake128 which has variable output length at this moment instead,
// since keccak256 hasher which has fixed output length (32 bytes)
use ecc::finite_field_arithmetic::bigint::BigInt;
use ecc::finite_field_arithmetic::traits::weierstrass_field::PrimeField;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};
use std::cmp::{Eq, PartialEq};
use std::marker::PhantomData;

const WORD_SIZE: usize = 64;
const NUM_LIMBS: usize = 4;
const STATE_SIZE: usize = NUM_LIMBS * (WORD_SIZE / 8);
type Word = u64;

#[derive(Clone, Debug)]
pub struct Shake128Transcript<F: PrimeField<NUM_LIMBS>> {
    pub round: Word,
    pub state: [u8; STATE_SIZE],
    hasher: Shake128,
    _p: PhantomData<F>,
}

impl<F: PrimeField<NUM_LIMBS>> Shake128Transcript<F> {
    fn hash(hasher: Shake128, input_bytes: &[u8]) -> [u8; STATE_SIZE] {
        let mut output_bytes = [0 as u8; STATE_SIZE];
        let mut local_hasher = hasher.clone();
        local_hasher.update(input_bytes);
        let mut reader = local_hasher.finalize_xof();
        reader.read(&mut output_bytes);
        output_bytes
    }

    // new a instance, initialize the state with a given input label
    // the hasher is still a empty one
    pub fn new(label: &'static [u8]) -> Self {
        let shake128_hasher = Shake128::default();
        let new_state = Self::hash(shake128_hasher.clone(), label);
        Self {
            round: 0 as Word,
            state: new_state,
            hasher: shake128_hasher,
            _p: Default::default(),
        }
    }

    // end of transcript session with a refreshed hasher
    pub fn squeeze(&mut self, label: &'static [u8]) -> F {
        // [round, state, label]
        let input_bytes: Vec<u8> = self
            .round
            .to_le_bytes()
            .into_iter()
            .chain(self.state.into_iter())
            .chain(label.to_vec())
            .collect();
        let new_state = Self::hash(self.hasher.clone(), input_bytes.as_slice());
        self.round = self.round + 1;
        self.state.copy_from_slice(&new_state);
        self.hasher = Shake128::default();

        // need to be reduced into scalar field
        F::from(BigInt::<NUM_LIMBS>::from(new_state.as_slice()))
    }

    // start of transcript session with a new hasher
    // absorb bytes, supporting any type of data
    pub fn absorb(&mut self, label: &'static [u8], input_bytes: &[u8]) {
        self.hasher.update(label);
        self.hasher.update(input_bytes);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ecc::finite_field_arithmetic::cycle_friendly::pallas::{fq::Fq, fr::Fr};
    use std::str::FromStr;
    type ScalarField = Fr<NUM_LIMBS>;
    type BaseField = Fq<NUM_LIMBS>;

    #[test]
    fn test_shake128() {
        let test_input = [b"10", b"20", b"30"];
        let mut all_bytes: Vec<u8> = Vec::new();

        // experiment object a
        let mut hasher_a = Shake128::default();
        for input in test_input {
            all_bytes.extend(input);
            hasher_a.update(input);
        }
        let mut reader_a = hasher_a.finalize_xof();
        let mut output_bytes_a = [0_u8; 4];
        reader_a.read(&mut output_bytes_a);

        // contradict object b
        let mut hasher_b = Shake128::default();
        hasher_b.update(all_bytes.as_slice());
        let mut reader_b = hasher_b.finalize_xof();
        let mut output_bytes_b = [0_u8; 4];
        reader_b.read(&mut output_bytes_b);

        assert_eq!(output_bytes_a, output_bytes_b);
    }

    #[test]
    fn test_transcript_single_session() {
        let mut all_bytes: Vec<u8> = Vec::new();

        // experiment object a
        let mut transcript_a = Shake128Transcript::<ScalarField>::new(b"TestInstance");
        // all_bytes.extend(b"TestInstance");
        let (init_round_bytes, init_state_bytes) =
            (transcript_a.round.to_le_bytes(), transcript_a.state);

        let test_input = ["10", "20", "30"];
        let test_label = [b"l1", b"l2", b"l3"];
        // start a transcript session
        for i in 0..test_input.len() {
            let scalar = ScalarField::from_str(test_input[i]).unwrap();
            // println!("{} scalar: {:?}", i, scalar);
            let scalar_bytes: [u8; STATE_SIZE] = scalar.to_bytes().try_into().unwrap();
            // println!("{} scalar bytes: {:?}", i, scalar_bytes);
            transcript_a.absorb(test_label[i], &scalar_bytes);
            all_bytes.extend(test_label[i]);
            all_bytes.extend(scalar_bytes);
        }
        // finish the transcript session
        let scalar_a = transcript_a.squeeze(b"output");

        // 2 * 3 + 32 * 3 + 6 + 8 + 32 = 148
        // contradict object b
        all_bytes.extend(init_round_bytes);
        all_bytes.extend(init_state_bytes);
        all_bytes.extend(b"output");
        assert_eq!(all_bytes.len(), 148);

        let mut hasher_b = Shake128::default();
        hasher_b.update(all_bytes.as_slice());
        let mut reader_b = hasher_b.finalize_xof();
        let mut output_bytes_b = [0 as u8; STATE_SIZE];
        reader_b.read(&mut output_bytes_b);
        let scalar_b = ScalarField::from(BigInt::<NUM_LIMBS>::from(output_bytes_b.as_slice()));

        // println!("---- output bytes: {:?}", output_bytes_b);

        assert_eq!(scalar_a, scalar_b);
    }

    #[test]
    fn test_multiple_session() {
        //////////////////// experiment object a
        let mut transcript_a = Shake128Transcript::<ScalarField>::new(b"TestInstance");
        // session 1
        let scalar_0_a = transcript_a.squeeze(b"init_state");
        transcript_a.absorb(b"label1", "123456789".as_bytes());
        transcript_a.absorb(b"label2", "987654321".as_bytes());
        let scalar_1_a = transcript_a.squeeze(b"session_1_state");

        // session 2
        transcript_a.absorb(b"label1", "123456789".as_bytes());
        transcript_a.absorb(b"label2", "987654321".as_bytes());
        let scalar_2_a = transcript_a.squeeze(b"session_2_state");

        //////////////////// contradiction object b
        // session 1
        let mut transcript_b = Shake128Transcript::<ScalarField>::new(b"TestInstance");
        let scalar_0_b = transcript_b.squeeze(b"init_state");
        transcript_b.absorb(b"label1", "123456789".as_bytes());
        transcript_b.absorb(b"label2", "987654321".as_bytes());
        let scalar_1_b = transcript_b.squeeze(b"session_1_state");

        // session 1
        transcript_b.absorb(b"label1", "123456789".as_bytes());
        transcript_b.absorb(b"label2", "987654321".as_bytes());
        let scalar_2_b = transcript_b.squeeze(b"session_2_state");
        assert_eq!(scalar_0_a, scalar_0_b);
        assert_eq!(scalar_1_a, scalar_1_b);
        assert_eq!(scalar_2_a, scalar_2_b);
    }
}
