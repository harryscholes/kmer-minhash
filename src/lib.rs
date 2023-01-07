use std::{
    collections::{hash_map::DefaultHasher, BinaryHeap},
    fmt::{self, Display, Formatter},
    hash::{Hash, Hasher},
};

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct Kmers<'a> {
    seq: &'a [u8],
    k: usize,
}

impl<'a> Kmers<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Kmers<'a> {
        Kmers { seq, k }
    }

    pub fn from_str(s: &'a str, k: usize) -> Kmers<'a> {
        Kmers {
            seq: s.as_bytes(),
            k,
        }
    }
}

impl<'a> IntoIterator for Kmers<'a> {
    type Item = &'a [u8];
    type IntoIter = KmersIntoIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        KmersIntoIter::new(self)
    }
}

#[derive(Clone, Copy)]
pub struct KmersIntoIter<'a> {
    kmers: Kmers<'a>,
    index: usize,
}

impl<'a> KmersIntoIter<'a> {
    pub fn new(kmers: Kmers<'a>) -> KmersIntoIter<'a> {
        KmersIntoIter { kmers, index: 0 }
    }

    pub fn new_from_index(kmers: Kmers<'a>, index: usize) -> KmersIntoIter<'a> {
        KmersIntoIter { kmers, index }
    }
}

impl<'a> Iterator for KmersIntoIter<'a> {
    type Item = &'a [u8];

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

pub trait MinHash<T>
where
    T: Hash,
{
    fn min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError>;
}

impl<'a> MinHash<&'a [u8]> for KmersIntoIter<'a> {
    fn min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError> {
        self.heap_min_hash(n)
    }
}

#[allow(dead_code)]
impl<'a> KmersIntoIter<'a> {
    fn heap_min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError> {
        let mut heap = BinaryHeap::with_capacity(n);

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
            Err(MinHashError::NotEnoughHashes { n: heap.len() })
        } else {
            Ok(heap.into_sorted_vec())
        }
    }
}

impl<'a> MinHash<&'a [u8]> for Kmers<'a> {
    fn min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError> {
        self.into_iter().min_hash(n)
    }
}

#[derive(Debug, PartialEq)]
pub enum MinHashError {
    NotEnoughHashes { n: usize },
}

impl std::error::Error for MinHashError {}

impl Display for MinHashError {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            MinHashError::NotEnoughHashes { n } => {
                write!(f, "{}", format!("Only {} hashes were produced", n))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::str;

    use bio::io::fasta;
    use rayon::prelude::*;
    use tokio::sync::mpsc;

    #[test]
    fn new() {
        let kmer = Kmers::new(&"abc".as_bytes(), 2);
        assert_eq!(
            kmer,
            Kmers {
                seq: &"abc".as_bytes(),
                k: 2
            }
        );
    }

    #[test]
    fn from_str() {
        let kmer = Kmers::from_str("abc", 2);
        assert_eq!(
            kmer,
            Kmers {
                seq: &"abc".as_bytes(),
                k: 2
            }
        );
    }

    #[test]
    fn iteration() {
        let kmers = Kmers::from_str("abcd", 2);
        let mut iter = kmers.into_iter();
        assert_eq!(iter.next().unwrap(), "ab".as_bytes());
        assert_eq!(iter.next().unwrap(), "bc".as_bytes());
        assert_eq!(iter.next().unwrap(), "cd".as_bytes());
        assert!(iter.next().is_none());
    }

    #[test]
    fn collect_iter() {
        let kmers = Kmers::from_str("abcd", 2);
        assert_eq!(
            kmers.into_iter().collect::<Vec<&[u8]>>(),
            vec!["ab".as_bytes(), "bc".as_bytes(), "cd".as_bytes()]
        );
    }

    #[test]
    fn kmers_iterator() {
        let kmers = Kmers::from_str("abcd", 2)
            .into_iter()
            .map(|s| s)
            .collect::<Vec<&[u8]>>();
        assert_eq!(
            kmers,
            vec!["ab".as_bytes(), "bc".as_bytes(), "cd".as_bytes()]
        )
    }

    #[test]
    fn kmers_iterator_from_index() {
        let kmers = KmersIntoIter::new_from_index(Kmers::from_str("abcd", 2), 1)
            .map(|s| s)
            .collect::<Vec<&[u8]>>();
        assert_eq!(kmers, vec!["bc".as_bytes(), "cd".as_bytes()])
    }

    #[test]
    fn min_hash() {
        let hashes = Kmers::from_str("abc", 2).min_hash(2).unwrap();
        assert_eq!(hashes.len(), 2);
        let mut manual_hashes = vec![hash("ab"), hash("bc")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes);
    }

    #[test]
    fn heap_min_hash() {
        let hashes = Kmers::from_str("abc", 2)
            .into_iter()
            .heap_min_hash(2)
            .unwrap();
        assert_eq!(hashes.len(), 2);
        let mut manual_hashes = vec![hash("ab"), hash("bc")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes);
    }

    #[test]
    fn min_hash_not_enough_hashes() {
        let err = Kmers::from_str("abcd", 2).min_hash(4).unwrap_err();
        assert_eq!(err, MinHashError::NotEnoughHashes { n: 3 });
    }

    #[test]
    fn sparse_min_hash_k2() {
        let hashes = Kmers::from_str("abcd", 2).min_hash(2).unwrap();
        let mut manual_hashes = vec![hash("ab"), hash("bc"), hash("cd")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..2]);
    }

    #[test]
    fn sparse_min_hash_k3() {
        let hashes = Kmers::from_str("abcdef", 3).min_hash(3).unwrap();
        let mut manual_hashes = vec![hash("abc"), hash("bcd"), hash("cde"), hash("def")];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..3]);
    }

    #[test]
    fn sparse_min_hash_k4() {
        let hashes = Kmers::from_str("abcdefgh", 4).min_hash(4).unwrap();
        let mut manual_hashes = vec![
            hash("abcd"),
            hash("bcde"),
            hash("cdef"),
            hash("defg"),
            hash("efgh"),
        ];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..4]);
    }

    #[test]
    fn concurrent_min_hash_rayon() {
        let kmers = vec![Kmers::from_str("abcd", 2), Kmers::from_str("bcde", 2)];
        let hashes = kmers
            .into_par_iter()
            .map(|kmers| kmers.min_hash(2).unwrap())
            .collect::<Vec<_>>();

        let mut manual_hashes_1 = vec![hash("ab"), hash("bc"), hash("cd")];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![hash("bc"), hash("cd"), hash("de")];
        manual_hashes_2.sort();
        let manual_hashes = vec![manual_hashes_1[..2].to_vec(), manual_hashes_2[..2].to_vec()];

        assert_eq!(hashes, manual_hashes);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn concurrent_min_hash_tokio() {
        let kmers = vec![Kmers::from_str("abcd", 2), Kmers::from_str("bcde", 2)];
        let mut all_hashes = vec![];

        let (tx, mut rx) = mpsc::unbounded_channel();
        for k in kmers {
            let tx = tx.clone();
            tokio::spawn(async move {
                let hashes = k.min_hash(2).unwrap();
                tx.send(hashes).unwrap();
            });
        }
        drop(tx);

        while let Some(hashes) = rx.recv().await {
            all_hashes.push(hashes);
        }

        let mut manual_hashes_1 = vec![hash("ab"), hash("bc"), hash("cd")];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![hash("bc"), hash("cd"), hash("de")];
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
                let kmers = Kmers::new(record.seq(), 3);
                let hashes = kmers.min_hash(3).unwrap();
                tx.send(hashes).unwrap();
            });
        }
        drop(tx);

        while let Some(hashes) = rx.recv().await {
            all_hashes.push(hashes);
        }

        let mut manual_hashes = vec!["abcdefgh", "ijklmnop"]
            .iter()
            .map(|s| Kmers::from_str(s, 3).into_iter().min_hash(3).unwrap())
            .collect::<Vec<Vec<u64>>>();

        all_hashes.sort();
        manual_hashes.sort();
        assert_eq!(all_hashes, manual_hashes);
    }

    fn hash(s: &str) -> u64 {
        let b = s.as_bytes();
        let mut h = DefaultHasher::new();
        b.hash(&mut h);
        h.finish()
    }
}
