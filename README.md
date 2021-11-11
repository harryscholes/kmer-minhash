# kmer-minhash

Parallel $k$-mer min-wise hashing in Rust.

This package provides:
- `Kmers` type for representing the $k$-mers of some sequence
- `KmersIterator` for iterating over the $k$-mers
- `MinHash` trait for min-wise hashing $k$-mers

### Example

```rust
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use kmer_minhash::{Kmers, MinHash};

const KMER_SIZE: usize = 3;
const N_HASHES: usize = 2;

fn main() {
    let seq = "abcde";

    let kmers = Kmers::new(seq, KMER_SIZE);

    let min_hashes = kmers.into_iter().min_hash(N_HASHES).unwrap();

    // Check that the minhashes are correct:
    let mut manual_hashes = vec!["abc", "bcd", "cde"]
        .iter()
        .map(|kmer| {
            let mut hasher = DefaultHasher::new();
            kmer.hash(&mut hasher);
            hasher.finish()
        })
        .collect::<Vec<u64>>();

    manual_hashes.sort();

    assert_eq!(min_hashes, manual_hashes[..N_HASHES]);
}
```

See the tests in `src/lib.rs` for more examples of how to use this package with Rayon and Tokio.