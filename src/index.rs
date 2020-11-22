use serde_derive::{Serialize, Deserialize};

use hashbrown::HashMap;

use debruijn::dna_string::DnaString;
use debruijn::Kmer;
use debruijn::msp::Scanner;

use std::cmp;

#[derive(Serialize, Deserialize)]
pub struct Index<P: Kmer> {
    reference: DnaString,
    map: HashMap<P, Vec<Interval>>,
    k: usize,
}

impl<P: Kmer> Index<P> {
    pub fn new(reference: &[u8], k: usize) -> Self {
        let reference = DnaString::from_acgt_bytes(reference);
        let score_func = |pi: &P| pi.to_u64() as usize;

        let scanner = Scanner::new(&reference, score_func, k);
        let intervals = scanner.scan(); // TODO: add to hash map while scanning to avoid extra allocations
        let mut map = HashMap::with_capacity(1 << 16);

        // map each minimizer p-mer to intervals of k-mers in the reference
        for i in &intervals {
            let bucket = map.entry(i.minimizer).or_insert_with(|| Vec::with_capacity(16));
            // TODO: use minimizer location to prune some cases when searching
            bucket.push(Interval::new(i.start, i.len));
        }

        Self { reference, map, k }
    }

    pub fn reference<'a>(&'a self) -> &'a DnaString {
        &self.reference
    }

    pub fn get<'a>(&'a self, p: &P) -> Option<&'a Vec<Interval>> {
        self.map.get(p)
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    pub fn k(&self) -> usize {
        self.k
    }

    /// Min, max, and avg len of vectors, and min, max, and avg len of all intervals.
    pub fn stats(&self) -> (usize, usize, f64, usize, usize, f64) {
        let mut min = usize::MAX;
        let mut max = usize::MIN;
        let mut sum = 0usize;
        let mut interval_sum_min = usize::MAX;
        let mut interval_sum_max = usize::MIN;
        let mut interval_sum_sum = 0usize;

        for v in self.map.values() {
            let len = v.len();
            min = cmp::min(min, len);
            max = cmp::max(max, len);
            sum += len;

            let mut interval_sum = 0usize;

            for i in v {
                interval_sum += i.len as usize;
            }

            interval_sum_min = cmp::min(interval_sum_min, interval_sum);
            interval_sum_max = cmp::max(interval_sum_max, interval_sum);
            interval_sum_sum += interval_sum;
        }

        let avg = (sum as f64) / (self.map.len() as f64);
        let interval_sum_avg = (interval_sum_sum as f64) / (self.map.len() as f64);
        (min, max, avg, interval_sum_min, interval_sum_max, interval_sum_avg)
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Interval {
    pub start: u32,
    pub len: u16 // TODO: may be able to use u8 as len?
}

impl Interval {
    pub fn new(start: u32, len: u16) -> Self {
        Self { start: start, len: len }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta;
    use debruijn::kmer::*;

    const CHR1_PATH: &'static str = "../dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa";

    #[test]
    fn test_chr1() {
        test_chr_index::<Kmer12>(CHR1_PATH, 32);
        test_chr_index::<Kmer16>(CHR1_PATH, 40);
    }

    fn test_chr_index<P: Kmer>(path: &str, k: usize) {
        let reader = fasta::Reader::from_file(path).unwrap();
        let mut num_records = 0usize;

        for record in reader.records() {
            let record = record.unwrap();

            let index: Index<P> = Index::new(record.seq(), k);
            let len = index.len();
            let stats = index.stats();

            println!("Record: {}\nMap len: {}\nMin bucket len: {}\nMax bucket len: {}\nAvg bucket len: {}", num_records, len, stats.0, stats.1, stats.2);
            println!("Min interval len per minimizer: {}\nMax interval len per minimizer: {}\nAvg interval len per minimizer: {}", stats.3, stats.4, stats.5);

            num_records += 1;
        }
    }
}
