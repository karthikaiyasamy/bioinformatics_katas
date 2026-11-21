use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

#[derive(Default, Debug, Clone)]
pub struct FastqMetrics {
    pub total_reads: u64,
    pub total_bases: u64,
    pub gc_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub quality_sum: u64,
}

impl FastqMetrics {
    pub fn merge(&mut self, other: &FastqMetrics) {
        self.total_reads += other.total_reads;
        self.total_bases += other.total_bases;
        self.gc_bases += other.gc_bases;
        self.q20_bases += other.q20_bases;
        self.q30_bases += other.q30_bases;
        self.quality_sum += other.quality_sum;
    }
}

/// Analyzes a batch of sequence & quality byte slices.
#[inline(always)]
pub fn analyze_record(seq: &[u8], qual: &[u8], metrics: &mut FastqMetrics) {
    let len = seq.len().min(qual.len());
    metrics.total_reads += 1;
    metrics.total_bases += len as u64;

    for i in 0..len {
        let base = seq[i];
        if base == b'G' || base == b'C' || base == b'g' || base == b'c' {
            metrics.gc_bases += 1;
        }

        // Standard Illumina Phred+33 decoding
        let q_score = (qual[i].saturating_sub(33)) as u64;
        metrics.quality_sum += q_score;

        if q_score >= 20 {
            metrics.q20_bases += 1;
        }
        if q_score >= 30 {
            metrics.q30_bases += 1;
        }
    }
}

/// Zero-allocation sequential streaming using a reusable buffer.
pub fn process_sequential(path: &Path) -> std::io::Result<FastqMetrics> {
    let file = File::open(path)?;
    // 1MB read buffer to minimize system call overhead
    let mut reader = BufReader::with_capacity(1024 * 1024, file);

    let mut metrics = FastqMetrics::default();
    let mut header = Vec::with_capacity(256);
    let mut seq = Vec::with_capacity(256);
    let mut plus = Vec::with_capacity(256);
    let mut qual = Vec::with_capacity(256);

    loop {
        header.clear();
        if reader.read_until(b'\n', &mut header)? == 0 {
            break; // EOF
        }

        seq.clear();
        if reader.read_until(b'\n', &mut seq)? == 0 {
            break;
        }

        plus.clear();
        if reader.read_until(b'\n', &mut plus)? == 0 {
            break;
        }

        qual.clear();
        if reader.read_until(b'\n', &mut qual)? == 0 {
            break;
        }

        // Trim trailing newline / carriage return bytes
        let seq_trimmed = trim_newline(&seq);
        let qual_trimmed = trim_newline(&qual);

        analyze_record(seq_trimmed, qual_trimmed, &mut metrics);
    }

    Ok(metrics)
}

/// Parallel batch processing using Rayon thread pool.
pub fn process_parallel(path: &Path, chunk_size: usize) -> std::io::Result<FastqMetrics> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(2 * 1024 * 1024, file);

    let mut overall_metrics = FastqMetrics::default();
    let mut batch: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(chunk_size);

    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut plus = Vec::new();
    let mut qual = Vec::new();

    loop {
        header.clear();
        if reader.read_until(b'\n', &mut header)? == 0 {
            break; // EOF
        }

        seq.clear();
        if reader.read_until(b'\n', &mut seq)? == 0 {
            break;
        }

        plus.clear();
        if reader.read_until(b'\n', &mut plus)? == 0 {
            break;
        }

        qual.clear();
        if reader.read_until(b'\n', &mut qual)? == 0 {
            break;
        }

        let seq_clean = trim_newline(&seq).to_vec();
        let qual_clean = trim_newline(&qual).to_vec();
        batch.push((seq_clean, qual_clean));

        if batch.len() >= chunk_size {
            let chunk_metrics = process_batch_rayon(&batch);
            overall_metrics.merge(&chunk_metrics);
            batch.clear();
        }
    }

    if !batch.is_empty() {
        let chunk_metrics = process_batch_rayon(&batch);
        overall_metrics.merge(&chunk_metrics);
    }

    Ok(overall_metrics)
}

fn process_batch_rayon(batch: &[(Vec<u8>, Vec<u8>)]) -> FastqMetrics {
    batch
        .par_iter()
        .fold(FastqMetrics::default, |mut acc, (seq, qual)| {
            analyze_record(seq, qual, &mut acc);
            acc
        })
        .reduce(FastqMetrics::default, |mut a, b| {
            a.merge(&b);
            a
        })
}

#[inline(always)]
fn trim_newline(slice: &[u8]) -> &[u8] {
    let mut end = slice.len();
    while end > 0 && (slice[end - 1] == b'\n' || slice[end - 1] == b'\r') {
        end -= 1;
    }
    &slice[..end]
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: fastq_summary_rs <path_to_fastq> [--parallel]");
        std::process::exit(1);
    }

    let fastq_path = Path::new(&args[1]);
    let is_parallel = args.iter().any(|arg| arg == "--parallel");

    let mode_str = if is_parallel { "Parallel (Rayon)" } else { "Sequential (Zero-Allocation)" };
    println!("=== Rust FASTQ Quality & GC Analyzer ===");
    println!("File: {:?}", fastq_path);
    println!("Mode: {}", mode_str);

    let start = Instant::now();
    let metrics = if is_parallel {
        process_parallel(fastq_path, 25_000).expect("Failed to process FASTQ in parallel")
    } else {
        process_sequential(fastq_path).expect("Failed to process FASTQ sequentially")
    };
    let duration = start.elapsed();

    let mean_qual = if metrics.total_bases > 0 {
        metrics.quality_sum as f64 / metrics.total_bases as f64
    } else {
        0.0
    };

    let gc_pct = if metrics.total_bases > 0 {
        (metrics.gc_bases as f64 / metrics.total_bases as f64) * 100.0
    } else {
        0.0
    };

    let q20_pct = if metrics.total_bases > 0 {
        (metrics.q20_bases as f64 / metrics.total_bases as f64) * 100.0
    } else {
        0.0
    };

    let q30_pct = if metrics.total_bases > 0 {
        (metrics.q30_bases as f64 / metrics.total_bases as f64) * 100.0
    } else {
        0.0
    };

    println!("\n--- Results ---");
    println!("Total Reads:        {:>12}", metrics.total_reads);
    println!("Total Bases:        {:>12}", metrics.total_bases);
    println!("Mean Quality (Q):   {:>12.2}", mean_qual);
    println!("GC Percentage:      {:>11.2}%", gc_pct);
    println!("Q20 Percentage:     {:>11.2}%", q20_pct);
    println!("Q30 Percentage:     {:>11.2}%", q30_pct);
    println!("Elapsed Time:       {:>12.3?}", duration);

    let throughput_mbs = (metrics.total_bases as f64 / 1_000_000.0) / duration.as_secs_f64();
    println!("Throughput:         {:>11.2} Megabases/sec", throughput_mbs);
}
