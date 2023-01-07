# Benchmark

Download some FASTA file, such as the [human proteome](https://www.uniprot.org/uniprot/?query=proteome:UP000005640%20reviewed:yes).

```console
$ cargo build --release
    Finished release [optimized] target(s) in 0.07s

$ time cat human.fasta | ../target/release/benchmark > human.minhash
cat human.fasta  0.00s user 0.01s system 10% cpu 0.097 total
../target/release/benchmark > human.minhash  0.60s user 0.10s system 370% cpu 0.189 total

$ wc human.*
  219680  422556 13604162 human.fasta
   20184   20184 25284945 human.minhash
```
