// for the purpose of research we use Shake128 which has variable output length at this moment instead,
// since keccak256 hasher which has fixed output length (32 bytes)
use ecc::finite_field_arithmetic::traits::{Field, PrimeField};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};
use std::marker::PhantomData;

const STATE_SIZE: usize = 2;

#[derive(Clone, Debug)]
pub struct Shake128Transcript<F: Field<STATE_SIZE>> {
    round: u8,
    state: [u8; STATE_SIZE],
    hasher: Shake128,
    _p: PhantomData<F>,
}

impl<F: Field<STATE_SIZE>> Shake128Transcript<F> {
    fn hash(hasher: Shake128, input_bytes: &[u8]) -> [u8; STATE_SIZE] {
        let mut output_bytes = [0u8; STATE_SIZE];
        let mut local_hasher = hasher.clone();
        local_hasher.update(input_bytes);
        let mut reader = local_hasher.finalize_xof();
        reader.read(&mut output_bytes);
        output_bytes
    }

    // new a instance, initialize the state with a given input label
    // the hasher is still a empty one
    fn new(label: &'static [u8]) -> Self {
        let shake128_hasher = Shake128::default();
        let new_state = Self::hash(shake128_hasher.clone(), label);
        Self {
            round: 0u8,
            state: new_state,
            hasher: shake128_hasher,
            _p: Default::default(),
        }
    }

    // end of transcript session with a refreshed hasher
    fn squeeze(&mut self, label: &'static [u8]) -> F {
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

        F::from(new_state)
    }

    // start of transcript session with a new hasher
    fn absorb(&mut self, label: &'static [u8], scalar: &F) {
        self.hasher.update(label);
        let scalar_bytes: [u8; STATE_SIZE] = (*scalar).into();
        self.hasher.update(&scalar_bytes);
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

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
        let mut transcript_a = Shake128Transcript::<PrimeField>::new(b"TestInstance");
        // all_bytes.extend(b"TestInstance");
        let (init_round_bytes, init_state_bytes) =
            (transcript_a.round.to_le_bytes(), transcript_a.state);

        let test_input = ["10", "20", "30"];
        let test_label = [b"l1", b"l2", b"l3"];
        // start a transcript session
        for i in 0..test_input.len() {
            let scalar = PrimeField::from_str(test_input[i]).unwrap();
            let scalar_bytes: [u8; 2] = scalar.into();
            transcript_a.absorb(test_label[i], &scalar);
            all_bytes.extend(test_label[i]);
            all_bytes.extend(scalar_bytes);
        }
        // finish the transcript session
        let scalar_a = transcript_a.squeeze(b"output");

        // contradict object b
        all_bytes.extend(init_round_bytes);
        all_bytes.extend(init_state_bytes);
        all_bytes.extend(b"output");

        let mut hasher_b = Shake128::default();
        hasher_b.update(all_bytes.as_slice());
        let mut reader_b = hasher_b.finalize_xof();
        let mut output_bytes_b = [0_u8; STATE_SIZE];
        reader_b.read(&mut output_bytes_b);
        let scalar_b = PrimeField::from(output_bytes_b);

        assert_eq!(scalar_a, scalar_b);
    }
}
