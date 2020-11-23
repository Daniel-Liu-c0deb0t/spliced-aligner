
// TODO: match any kmer instead of sequentially calculating the MMP!
// using a kmer-based index means that short matches cannot be found, especially when there are
// splice junctions

pub fn max_map_prefix<K: Kmer, P: Kmer>(query: &DnaString, index: &Index<K, P>) -> Option<(u32, u32)> {
    let kmer = query.get_kmer::<K>(0);
    let minimizer: P = minimizer(query.prefix(index.k()));
    let intervals = index.get(minimizer);
    let reference = index.reference();
    let max_pos = None;
    let max_len = 0usize;

    for i in intervals {
        let mut interval_kmer = reference.get_kmer(i.start);
        let end = i.start + i.len + 1 - K::k();

        for start in i.start..end {
            if interval_kmer == kmer {
                let len = K::k() + lcp(reference.slice(start + K::k(), reference.len()), query);

                if len > max_len {
                    max_len = len;
                    max_pos = Some(start);
                }
            }

            if start < end - 1 {
                interval_kmer = interval_kmer.extend_right(reference.get(start + K::k()));
            }
        }
    }

    (max_pos?, max_len)
}

fn minimizer<V: Vmer, P: Kmer>(v: V) -> P {
    let mut min = v.get_kmer::<P>(0);
    let mut curr = min.clone();

    for i in 1..(v.len() + 1 - P::k()) {
        curr = curr.extend_right(v.get(i + P::k() - 1));
        min = cmp::min(min, curr);
    }

    min
}

/// Longest common prefix.
fn lcp<V1: Vmer, V2: Vmer>(a: V1, b: V2) -> usize {
    let min_len = cmp::min(a.len(), b.len());
    let mut idx = 0;
    let k = Kmer32::k();

    // match 32 bp at once
    while idx + k <= min_len {
        let kmer1 = a.get_kmer::<Kmer32>(idx);
        let kmer2 = b.get_kmer::<Kmer32>(idx);

        if(kmer1 != kmer2)
            break;

        idx += k;
    }

    // TODO: potentially replaced by a fast kmer LCP function?
    while idx < min_len {
        if a.get(idx) != b.get(idx) {
            break;
        }

        idx++;
    }

    idx
}
