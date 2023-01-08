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
        Self::at_index(kmers, 0)
    }

    pub fn at_index(kmers: Kmers<'a>, index: usize) -> KmersIterator<'a> {
        KmersIterator { kmers, index }
    }
}

impl<'a> Iterator for KmersIterator<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        let end = self.index + self.kmers.k;

        if end <= self.kmers.seq.len() {
            let kmer = &self.kmers.seq[self.index..end];
            self.index += 1;
            Some(kmer)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let l = 1 + self.kmers.seq.len() - self.index - self.kmers.k;
        (l, Some(l))
    }
}

impl<'a> std::iter::ExactSizeIterator for KmersIterator<'a> {}

impl<'a> std::iter::FusedIterator for KmersIterator<'a> {}

pub trait MinHash<T>
where
    T: Hash,
{
    fn min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError>;
}

impl<'a> MinHash<&'a [u8]> for Kmers<'a> {
    fn min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError> {
        self.heap_min_hash(n)
    }
}

impl<'a> Kmers<'a> {
    fn heap_min_hash(&self, n: usize) -> Result<Vec<u64>, MinHashError> {
        let mut heap = BinaryHeap::with_capacity(n);

        for kmer in self.into_iter() {
            let hash = hash(kmer);

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
        }

        if heap.len() < n {
            Err(MinHashError::NotEnoughHashes { n: heap.len() })
        } else {
            Ok(heap.into_sorted_vec())
        }
    }
}

fn hash(input: impl Hash) -> u64 {
    let mut hasher = DefaultHasher::new();
    input.hash(&mut hasher);
    hasher.finish()
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

    use bio::io::fasta;
    use rayon::prelude::*;
    use tokio::sync::mpsc;

    #[test]
    fn new() {
        let kmer = Kmers::new("abc".as_bytes(), 2);
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
    fn kmers_iterator_at_index() {
        let kmers = KmersIterator::at_index(Kmers::from_str("abcd", 2), 1)
            .map(|s| s)
            .collect::<Vec<&[u8]>>();
        assert_eq!(kmers, vec!["bc".as_bytes(), "cd".as_bytes()])
    }

    #[test]
    fn kmers_iterator_len_k2() {
        let mut iter = KmersIterator::new(Kmers::from_str("abcd", 2));
        for i in 3..=0 {
            assert_eq!(iter.len(), i);
            iter.next();
        }
    }

    #[test]
    fn kmers_iterator_len_k3() {
        let mut iter = KmersIterator::new(Kmers::from_str("abcde", 3));
        for i in 3..=0 {
            assert_eq!(iter.len(), i);
            iter.next();
        }
    }

    #[test]
    fn min_hash() {
        let kmers = Kmers::from_str("abc", 2);
        let hashes = kmers.min_hash(2).unwrap();
        assert_eq!(hashes.len(), 2);
        let mut manual_hashes = vec![hash("ab".as_bytes()), hash("bc".as_bytes())];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes);
    }

    #[test]
    fn heap_min_hash() {
        let hashes = Kmers::from_str("abc", 2).heap_min_hash(2).unwrap();
        assert_eq!(hashes.len(), 2);
        let mut manual_hashes = vec![hash("ab".as_bytes()), hash("bc".as_bytes())];
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
        let mut manual_hashes = vec![
            hash("ab".as_bytes()),
            hash("bc".as_bytes()),
            hash("cd".as_bytes()),
        ];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..2]);
    }

    #[test]
    fn sparse_min_hash_k3() {
        let hashes = Kmers::from_str("abcdef", 3).min_hash(3).unwrap();
        let mut manual_hashes = vec![
            hash("abc".as_bytes()),
            hash("bcd".as_bytes()),
            hash("cde".as_bytes()),
            hash("def".as_bytes()),
        ];
        manual_hashes.sort();
        assert_eq!(hashes, manual_hashes[..3]);
    }

    #[test]
    fn sparse_min_hash_k4() {
        let hashes = Kmers::from_str("abcdefgh", 4).min_hash(4).unwrap();
        let mut manual_hashes = vec![
            hash("abcd".as_bytes()),
            hash("bcde".as_bytes()),
            hash("cdef".as_bytes()),
            hash("defg".as_bytes()),
            hash("efgh".as_bytes()),
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

        let mut manual_hashes_1 = vec![
            hash("ab".as_bytes()),
            hash("bc".as_bytes()),
            hash("cd".as_bytes()),
        ];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![
            hash("bc".as_bytes()),
            hash("cd".as_bytes()),
            hash("de".as_bytes()),
        ];
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

        let mut manual_hashes_1 = vec![
            hash("ab".as_bytes()),
            hash("bc".as_bytes()),
            hash("cd".as_bytes()),
        ];
        manual_hashes_1.sort();
        let mut manual_hashes_2 = vec![
            hash("bc".as_bytes()),
            hash("cd".as_bytes()),
            hash("de".as_bytes()),
        ];
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
            .map(|s| Kmers::from_str(s, 3).min_hash(3).unwrap())
            .collect::<Vec<Vec<u64>>>();

        all_hashes.sort();
        manual_hashes.sort();
        assert_eq!(all_hashes, manual_hashes);
    }
}
