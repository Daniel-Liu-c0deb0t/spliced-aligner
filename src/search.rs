
// TODO: match any kmer instead of sequentially calculating the MMP!
// using a kmer-based index means that short matches cannot be found, especially when there are
// splice junctions

pub fn max_map_prefix<P: Kmer>(query: &DnaString, index: &Index<P>) -> Option<(u32, u32)> {
    // TODO: use Kmer instead of DnaStringSlice (Index includes Kmer type)
    let kmer = query.prefix(index.k());
    let minimizer = minimizer(kmer);
    let intervals = index.get(minimizer);
    let reference = index.reference();
    let max_pos = None;
    let max_len = 0usize;

    for i in intervals {
        for start in i.start..(i.start + i.len + 1 - index.k()) {
            // TODO: if there is a kmer function for fast LCP, then change this
            let slice = reference.slice(start, start + index.k());

            if slice == kmer {
                let len = index.k() + lcp(reference.slice(start + index.k(), reference.len()), query);

                if len > max_len {
                    max_len = len;
                    max_pos = Some(start);
                }
            }
        }
    }

    (max_pos?, max_len)
}

fn minimizer<P: Kmer, V: Vmer>(v: V) -> P {

}
