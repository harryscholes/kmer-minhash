use std::io::{self, Write};

use bio::io::fasta;
use itertools::Itertools;
use tokio::sync::mpsc;

use kmer_minhash::{Kmers, MinHash, MinHashError};

const KMER_SIZE: usize = 3;
const N_HASHES: usize = 64;

#[tokio::main]
async fn main() {
    let reader = fasta::Reader::new(io::stdin());

    let (tx, mut rx) = mpsc::unbounded_channel();
    for record_result in reader.records() {
        let tx = tx.clone();
        tokio::spawn(async move {
            let record = record_result.unwrap();
            let kmers = Kmers::new(record.seq(), KMER_SIZE);
            match kmers.into_iter().min_hash(N_HASHES) {
                Ok(hashes) => tx.send(hashes).unwrap(),
                Err(e) => match e {
                    MinHashError::NotEnoughHashes { .. } => {
                        eprintln!("{} for sequence id {}", e, record.id())
                    }
                },
            }
        });
    }
    drop(tx);

    let stdout = std::io::stdout();
    let mut lock = stdout.lock();
    while let Some(hashes) = rx.recv().await {
        writeln!(lock, "{}", hashes.iter().join(",")).unwrap();
    }
}
