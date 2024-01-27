type Word = u64;
const WORD_SIZE: usize = 64;

// little-endian
pub fn bytes_to_word(bytes: &[u8]) -> Word {
    assert!(bytes.len() < WORD_SIZE);
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
