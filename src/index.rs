use serde_derive::{Serialize, Deserialize};

use hashbrown::HashMap;

use debruijn::dna_string::DnaString;
use debruijn::Kmer;
use debruijn::msp::Scanner;

#[derive(Serialize, Deserialize)]
pub struct Index<P: Kmer> {
    reference: DnaString,
    map: HashMap<P, Vec<Interval>>
}

impl<P: Kmer> Index<P> {
    pub fn new(reference: &[u8], k: usize) -> Self {
        let sequence = DnaString::from_acgt_bytes(reference);
        let score_func = |pi: &P| pi.to_u64() as usize;

        let scanner = Scanner::new(&sequence, score_func, k);
        let intervals = scanner.scan(); // TODO: add to hash map while scanning to avoid extra allocations
        let mut map = HashMap::with_capacity(1 << 20);

        // map each minimizer p-mer to intervals of k-mers in the reference
        for i in &intervals {
            let bucket = map.entry(i.minimizer).or_insert_with(|| Vec::with_capacity(16));
            // TODO: use minimizer location to prune some cases when searching
            bucket.push(Interval::new(i.start, i.len));
        }

        Self {
            reference: sequence,
            map: map,
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Interval {
    start: u32,
    len: u16
}

impl Interval {
    fn new(start: u32, len: u16) -> Self {
        Self { start: start, len: len }
    }
}
