use std::{
    io::{self, Write},
    str,
};

use bio::io::fasta;
use itertools::Itertools;
use tokio::sync::mpsc;

use kmer_minhash::{Kmers, MinHash};

#[tokio::main]
async fn main() {
    let reader = fasta::Reader::new(io::stdin());

    let (tx, mut rx) = mpsc::unbounded_channel();
    for record_result in reader.records() {
        let tx = tx.clone();
        tokio::spawn(async move {
            let record = record_result.unwrap();
            let seq = str::from_utf8(record.seq()).unwrap();
            let k = Kmers::new(seq, 8);
            let hashes = k.into_iter().min_hash(64).unwrap();
            tx.send(hashes).unwrap();
        });
    }
    drop(tx);

    let stdout = std::io::stdout();
    let mut lock = stdout.lock();
    while let Some(hashes) = rx.recv().await {
        writeln!(lock, "{}", hashes.iter().join(",")).unwrap();
    }
}
