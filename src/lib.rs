use std::collections::hash_map::DefaultHasher;
use std::collections::BinaryHeap;
use std::hash::{Hash, Hasher};

// TODO implement with `seq` as `&[u8]`
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct Kmers<'a> {
    seq: &'a str,
    k: usize,
}

impl<'a> Kmers<'a> {
    pub fn new(seq: &'a str, k: usize) -> Kmers<'a> {
        Kmers { seq, k }
    }
}

impl<'a> IntoIterator for Kmers<'a> {
    type Item = &'a str;
    type IntoIter = KmersIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        KmersIterator::new(self)
    }
}

pub struct KmersIterator<'a> {
    kmers: Kmers<'a>,
    index: usize,
}

impl<'a> KmersIterator<'a> {
    pub fn new(kmers: Kmers<'a>) -> KmersIterator<'a> {
        KmersIterator { kmers, index: 0 }
    }

    pub fn new_from_index(kmers: Kmers<'a>, index: usize) -> KmersIterator<'a> {
        KmersIterator { kmers, index }
    }
}

impl<'a> Iterator for KmersIterator<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        let (start, end) = (self.index, self.index + self.kmers.k);
        if end > self.kmers.seq.len() {
            None
        } else {
            self.index += 1;
            Some(&self.kmers.seq[start..end])
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum MinHashError {
    NotEnoughHashes,
}

use std::fmt::{self, Display, Formatter};

impl Display for MinHashError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Not enough hashes were produced")
    }
}

impl std::error::Error for MinHashError {}

pub trait MinHash<T>: Iterator<Item = T>
where
    T: Hash,
{
    fn min_hash(&mut self, n: usize) -> Result<Vec<u64>, MinHashError> {
        let mut heap = BinaryHeap::new();

        self.for_each(|s| {
            let mut hasher = DefaultHasher::new();
            s.hash(&mut hasher);
            let hash = hasher.finish();

            match heap.peek() {
                None => heap.push(hash),
                Some(last) => {
                    if heap.len() < n {
                        heap.push(hash);
                    } else if hash < *last {
                        heap.pop();
                        heap.push(hash);
                    }
                }
            }
        });

        if heap.len() < n {
            Err(MinHashError::NotEnoughHashes {})
        } else {
            Ok(heap.into_sorted_vec())
        }
    }
}

impl<'a> MinHash<&'a str> for KmersIterator<'a> {}

#[cfg(test)]
mod tests {
    use super::*;

    use std::str;

    use bio::io::fasta;
    use rayon::prelude::*;
    use tokio::sync::mpsc;

    #[test]
    fn constructor() {
        let kmer = Kmers::new("abc", 2);
        assert_eq!(kmer, Kmers { seq: "abc", k: 2 });
    }

    #[test]
    fn iteration() {
        let kmer = Kmers::new("abcd", 2);
        let mut iter = kmer.into_iter();
        assert_eq!(iter.next().unwrap(), "ab");
        assert_eq!(iter.next().unwrap(), "bc");
        assert_eq!(iter.next().unwrap(), "cd");
        assert!(iter.next().is_none());
    }

    #[test]
    fn iterator() {
        let kmers = Kmers::new("abcd", 2)
            .into_iter()
            .map(|s| s)
            .collect::<Vec<&str>>();
        assert_eq!(kmers, vec!["ab", "bc", "cd"])
    }

    #[test]
    fn iterator_from_index() {
        let kmers = KmersIterator::new_from_index(Kmers::new("abcd", 2), 1)
            .map(|s| s)
            .collect::<Vec<&str>>();
        assert_eq!(kmers, vec!["bc", "cd"])
    }

    #[test]
    fn min_hash() {
        let hashes = Kmers::new("abc", 2).into_iter().min_hash(2).unwrap();
        assert_eq!(hashes.len(), 2);
        let mut manual_hashes = vec![hash_str("ab"), hash_str("bc")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes);
    }

    #[test]
    fn not_enough_hashes() {
        let err = Kmers::new("abcd", 2).into_iter().min_hash(4).unwrap_err();
        assert_eq!(err, MinHashError::NotEnoughHashes {});
    }

    #[test]
    fn sparse_min_hash_k2() {
        let hashes = Kmers::new("abcd", 2).into_iter().min_hash(2).unwrap();
        let mut manual_hashes = vec![hash_str("ab"), hash_str("bc"), hash_str("cd")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..2]);
    }

    #[test]
    fn sparse_min_hash_k3() {
        let hashes = Kmers::new("abcdef", 3).into_iter().min_hash(3).unwrap();
        let mut manual_hashes = vec![
            hash_str("abc"),
            hash_str("bcd"),
            hash_str("cde"),
            hash_str("def"),
        ];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..3]);
    }

    #[test]
    fn sparse_min_hash_k4() {
        let hashes = Kmers::new("abcdefgh", 4).into_iter().min_hash(4).unwrap();
        let mut manual_hashes = vec![
            hash_str("abcd"),
            hash_str("bcde"),
            hash_str("cdef"),
            hash_str("defg"),
            hash_str("efgh"),
        ];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..4]);
    }

    #[test]
    fn concurrent_min_hash_rayon() {
        let kmers = vec![Kmers::new("abcd", 2), Kmers::new("bcde", 2)];
        let hashes = kmers
            .par_iter()
            .map(|k| k.into_iter().min_hash(2).unwrap())
            .collect::<Vec<_>>();

        let mut manual_hashes_1 = vec![hash_str("ab"), hash_str("bc"), hash_str("cd")];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![hash_str("bc"), hash_str("cd"), hash_str("de")];
        manual_hashes_2.sort();
        let manual_hashes = vec![manual_hashes_1[..2].to_vec(), manual_hashes_2[..2].to_vec()];

        assert_eq!(hashes, manual_hashes);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn concurrent_min_hash_tokio() {
        let kmers = vec![Kmers::new("abcd", 2), Kmers::new("bcde", 2)];
        let mut all_hashes = vec![];

        let (tx, mut rx) = mpsc::unbounded_channel();
        for k in kmers {
            let tx = tx.clone();
            tokio::spawn(async move {
                let hashes = k.into_iter().min_hash(2).unwrap();
                tx.send(hashes).unwrap();
            });
        }
        drop(tx);

        while let Some(hashes) = rx.recv().await {
            all_hashes.push(hashes);
        }

        let mut manual_hashes_1 = vec![hash_str("ab"), hash_str("bc"), hash_str("cd")];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![hash_str("bc"), hash_str("cd"), hash_str("de")];
        manual_hashes_2.sort();
        let mut manual_hashes = vec![manual_hashes_1[..2].to_vec(), manual_hashes_2[..2].to_vec()];

        all_hashes.sort();
        manual_hashes.sort();
        assert_eq!(all_hashes, manual_hashes);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn concurrent_min_hash_tokio_fasta() {
        let fasta_file: &[u8] = b"> seq_1
abcd
efgh
> seq_2
ijkl
mnop
";
        let reader = fasta::Reader::new(fasta_file);
        let mut all_hashes = vec![];

        let (tx, mut rx) = mpsc::unbounded_channel();
        for record_result in reader.records() {
            let tx = tx.clone();
            tokio::spawn(async move {
                let record = record_result.unwrap();
                let seq = str::from_utf8(record.seq()).unwrap();
                let k = Kmers::new(seq, 3);
                let hashes = k.into_iter().min_hash(3).unwrap();
                tx.send(hashes).unwrap();
            });
        }
        drop(tx);

        while let Some(hashes) = rx.recv().await {
            all_hashes.push(hashes);
        }

        let mut manual_hashes = vec!["abcdefgh", "ijklmnop"]
            .iter()
            .map(|s| Kmers::new(s, 3).into_iter().min_hash(3).unwrap())
            .collect::<Vec<Vec<u64>>>();

        all_hashes.sort();
        manual_hashes.sort();
        assert_eq!(all_hashes, manual_hashes);
    }

    fn hash_str(s: &str) -> u64 {
        let mut h = DefaultHasher::new();
        s.hash(&mut h);
        h.finish()
    }
}
