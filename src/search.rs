use debruijn::dna_string::DnaStringSlice;
use debruijn::{Kmer, Vmer, Mer};
use debruijn::kmer::Kmer32;

use std::cmp;

use super::index::Index;

// TODO: match any kmer instead of sequentially calculating the MMP!
// using a kmer-based index means that short matches cannot be found, especially when there are
// splice junctions

/// Returns an optional (start, length) tuple for the longest mappable location.
pub fn max_map_prefix<K: Kmer, P: Kmer>(index: &Index<K, P>, query: DnaStringSlice) -> Option<(usize, usize)> {
    let kmer = query.get_kmer::<K>(0);
    let minimizer: P = minimizer(&query.slice(0, K::k()));
    let query_end = query.slice(K::k(), query.len());
    let intervals = index.get(minimizer)?;
    let reference = index.reference();
    let mut max_pos = None;
    let mut max_len = 0usize;

    for i in intervals {
        let interval_start = i.start as usize;
        let interval_len = i.len as usize;
        let mut interval_kmer = reference.get_kmer::<K>(interval_start).extend_left(0);
        let end = interval_start + interval_len + 1 - K::k();

        for start in interval_start..end {
            interval_kmer = interval_kmer.extend_right(reference.get(start + K::k() - 1));

            if interval_kmer == kmer {
                // if kmer match, then try extend
                let lcp = lcp(&reference.slice(start + K::k(), reference.len()), &query_end);
                let len = K::k() + lcp;

                if len > max_len {
                    max_len = len;
                    max_pos = Some(start);
                }
            }
        }
    }

    Some((max_pos?, max_len))
}

fn minimizer<V: Vmer, P: Kmer>(v: &V) -> P {
    let mut min = v.get_kmer::<P>(0);
    let mut curr = min.clone();

    for i in 1..(v.len() + 1 - P::k()) {
        curr = curr.extend_right(v.get(i + P::k() - 1));
        min = cmp::min(min, curr);
    }

    min
}

/// Longest common prefix.
fn lcp<V1: Vmer, V2: Vmer>(a: &V1, b: &V2) -> usize {
    let min_len = cmp::min(a.len(), b.len());
    let mut idx = 0;
    let k = Kmer32::k();

    // match 32 bp at once
    while idx + k <= min_len {
        let kmer1 = a.get_kmer::<Kmer32>(idx);
        let kmer2 = b.get_kmer::<Kmer32>(idx);

        if kmer1 != kmer2 {
            break;
        }

        idx += k;
    }

    // TODO: potentially replaced by a fast kmer LCP function?
    while idx < min_len {
        if a.get(idx) != b.get(idx) {
            break;
        }

        idx += 1;
    }

    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::*;
    use bio::io::fasta;
    use debruijn::kmer::*;
    use debruijn::dna_string::*;

    #[test]
    fn test_chr1() {
        test_chr_mmp::<Kmer32, Kmer12>(crate::tests::CHR1_PATH,
                                       b"TTGCGAACACCAAGCTCAACAATGAGCCCTGGAAAATTTCTGGAATGGATTATTAAACAG",
                                       Some((334 * 60, 60))); // fasta file line 336

        test_chr_mmp::<Kmer32, Kmer12>(crate::tests::CHR1_PATH,
                                       b"TTGCGAACACCAAGCTCAACAATGAGCCCTGGAAAATTTCTGGAATGGATCCCCCCCCCC",
                                       Some((334 * 60, 50))); // mismatching tail
    }

    fn test_chr_mmp<K: Kmer, P: Kmer>(path: &str, query: &[u8], correct: Option<(usize, usize)>) {
        let query = DnaString::from_acgt_bytes(query);
        let reader = fasta::Reader::from_file(path).unwrap();

        for record in reader.records() {
            let record = record.unwrap();

            let index: Index<K, P> = Index::new(record.seq());

            println!("Index building done!");

            let query_slice = query.slice(0, query.len());
            let res = max_map_prefix::<K, P>(&index, query_slice);

            println!("Query: {:?}", query);
            println!("MMP location: {:?}", res);

            assert_eq!(correct, res);
        }
    }
}
