//! Biometal PCR Duplicate Detection Tool
//!
//! Fast PCR duplicate detection using biometal k-mer primitives
//! Uses minimizer-based sequence fingerprinting to replace clumpify PCR deduplication

use anyhow::Result;
use biometal::operations::kmer::{extract_minimizers_fast, Minimizer};
use biometal::{FastqStream, FastqWriter, FastqRecord};
use clap::{Arg, Command};
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashSet;
use std::path::PathBuf;

#[derive(Debug, Serialize, Deserialize)]
struct PcrDedupStats {
    total_reads: u64,
    pcr_duplicates_found: u64,
    unique_reads_kept: u64,
    duplicate_clusters: u64,
    average_cluster_size: f64,
    similarity_threshold: f64,
    kmer_size: usize,
    window_size: usize,
    processing_time_seconds: f64,
}

fn main() -> Result<()> {
    let matches = Command::new("biometal-pcr-dedup")
        .version("0.1.0")
        .about("Fast PCR duplicate detection using biometal k-mer primitives")
        .author("Megan Johnson")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FASTQ")
                .help("Input FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FASTQ")
                .help("Output deduplicated FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output PCR deduplication statistics (JSON)")
                .default_value("pcr_dedup_stats.json"),
        )
        .arg(
            Arg::new("threshold")
                .short('t')
                .long("threshold")
                .value_name("SIMILARITY")
                .help("Jaccard similarity threshold for PCR duplicates (0.0-1.0)")
                .default_value("0.8"),
        )
        .arg(
            Arg::new("kmer_size")
                .short('k')
                .long("kmer-size")
                .value_name("SIZE")
                .help("K-mer size for minimizer extraction")
                .default_value("21"),
        )
        .arg(
            Arg::new("window_size")
                .short('w')
                .long("window-size")
                .value_name("SIZE")
                .help("Window size for minimizer extraction")
                .default_value("11"),
        )
        .arg(
            Arg::new("keep_best")
                .long("keep-best-quality")
                .help("Keep highest quality read from each duplicate cluster")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    // Parse arguments
    let input_path: PathBuf = matches.get_one::<String>("input").unwrap().into();
    let output_path: PathBuf = matches.get_one::<String>("output").unwrap().into();
    let stats_path: PathBuf = matches.get_one::<String>("stats").unwrap().into();
    let threshold: f64 = matches.get_one::<String>("threshold").unwrap().parse()?;
    let kmer_size: usize = matches.get_one::<String>("kmer_size").unwrap().parse()?;
    let window_size: usize = matches.get_one::<String>("window_size").unwrap().parse()?;
    let keep_best_quality = matches.get_flag("keep_best");

    // Validate parameters
    if threshold < 0.0 || threshold > 1.0 {
        return Err(anyhow::anyhow!("Similarity threshold must be between 0.0 and 1.0"));
    }

    println!("🧬 Biometal PCR Duplicate Detection");
    println!("==================================");
    println!("Input: {}", input_path.display());
    println!("Output: {}", output_path.display());
    println!("Similarity threshold: {:.2}", threshold);
    println!("K-mer size: {}, Window size: {}", kmer_size, window_size);
    println!("Quality selection: {}", if keep_best_quality { "Best quality" } else { "First occurrence" });

    let start_time = std::time::Instant::now();

    // Process PCR duplicates
    let stats = process_pcr_duplicates(
        &input_path,
        &output_path,
        threshold,
        kmer_size,
        window_size,
        keep_best_quality,
    )?;

    let processing_time = start_time.elapsed().as_secs_f64();
    let final_stats = PcrDedupStats {
        processing_time_seconds: processing_time,
        ..stats
    };

    // Write statistics
    let stats_json = serde_json::to_string_pretty(&final_stats)?;
    std::fs::write(&stats_path, stats_json)?;

    println!("\n✅ PCR Deduplication Complete");
    println!("Total reads processed: {}", final_stats.total_reads);
    println!("PCR duplicates found: {}", final_stats.pcr_duplicates_found);
    println!("Unique reads kept: {}", final_stats.unique_reads_kept);
    println!("Duplicate clusters: {}", final_stats.duplicate_clusters);
    println!("Average cluster size: {:.2}", final_stats.average_cluster_size);
    println!("Processing time: {:.2}s", processing_time);
    println!("Statistics written to: {}", stats_path.display());

    Ok(())
}

fn process_pcr_duplicates(
    input_path: &PathBuf,
    output_path: &PathBuf,
    threshold: f64,
    kmer_size: usize,
    window_size: usize,
    keep_best_quality: bool,
) -> Result<PcrDedupStats> {
    // Step 1: Read all records and extract minimizer signatures
    let fastq_stream = FastqStream::from_path(input_path)?;
    let mut records_with_signatures: Vec<(FastqRecord, Vec<Minimizer>)> = Vec::new();
    let mut total_reads = 0;
    let mut signature_errors = 0;

    println!("📊 Reading FASTQ and extracting minimizer signatures...");

    for record_result in fastq_stream {
        let record = record_result?;
        total_reads += 1;

        match extract_minimizers_fast(&record.sequence, kmer_size, window_size) {
            Ok(minimizers) => {
                records_with_signatures.push((record, minimizers));
            }
            Err(_) => {
                signature_errors += 1;
                // Keep records that can't be processed with empty signature
                records_with_signatures.push((record, Vec::new()));
            }
        }

        if total_reads % 10000 == 0 {
            println!("   - Processed {} reads...", total_reads);
        }
    }

    println!("   - Total reads: {}", total_reads);
    println!("   - Successful signatures: {}", records_with_signatures.len() - signature_errors);
    println!("   - Signature errors: {}", signature_errors);

    // Step 2: Calculate pairwise similarities and cluster duplicates
    println!("🔍 Calculating sequence similarities...");

    let duplicate_clusters = find_pcr_duplicates(&records_with_signatures, threshold)?;

    println!("   - Found {} duplicate clusters", duplicate_clusters.len());

    // Step 3: Create duplicate index mapping
    let mut duplicate_indices: HashSet<usize> = HashSet::new();
    let mut cluster_sizes: Vec<usize> = Vec::new();

    for cluster in &duplicate_clusters {
        cluster_sizes.push(cluster.len());

        if cluster.len() > 1 {
            // Select representative (first or best quality)
            let representative = if keep_best_quality {
                select_best_quality_read(cluster, &records_with_signatures)?
            } else {
                cluster[0] // First occurrence
            };

            // Mark all others as duplicates
            for &idx in cluster {
                if idx != representative {
                    duplicate_indices.insert(idx);
                }
            }
        }
    }

    println!("   - PCR duplicates to remove: {}", duplicate_indices.len());

    // Step 4: Write filtered output
    println!("📝 Writing deduplicated output...");

    let mut writer = FastqWriter::create(output_path)?;
    let mut unique_reads_kept = 0;

    for (i, (record, _)) in records_with_signatures.iter().enumerate() {
        if !duplicate_indices.contains(&i) {
            writer.write_record(record)?;
            unique_reads_kept += 1;
        }
    }

    let average_cluster_size = if duplicate_clusters.is_empty() {
        0.0
    } else {
        cluster_sizes.iter().sum::<usize>() as f64 / duplicate_clusters.len() as f64
    };

    Ok(PcrDedupStats {
        total_reads,
        pcr_duplicates_found: duplicate_indices.len() as u64,
        unique_reads_kept,
        duplicate_clusters: duplicate_clusters.len() as u64,
        average_cluster_size,
        similarity_threshold: threshold,
        kmer_size,
        window_size,
        processing_time_seconds: 0.0, // Will be set by caller
    })
}

fn find_pcr_duplicates(
    records_with_signatures: &[(FastqRecord, Vec<Minimizer>)],
    threshold: f64,
) -> Result<Vec<Vec<usize>>> {
    let n_records = records_with_signatures.len();
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut assigned: HashSet<usize> = HashSet::new();

    // Process in batches to manage memory for large datasets
    let batch_size = std::cmp::min(1000, n_records);

    for batch_start in (0..n_records).step_by(batch_size) {
        let batch_end = std::cmp::min(batch_start + batch_size, n_records);

        if batch_start > 0 && batch_start % 10000 == 0 {
            println!("   - Similarity analysis: {}/{} records processed", batch_start, n_records);
        }

        for i in batch_start..batch_end {
            if assigned.contains(&i) {
                continue; // Already in a cluster
            }

            let mut current_cluster = vec![i];
            assigned.insert(i);

            // Compare with all subsequent records
            for j in (i + 1)..n_records {
                if assigned.contains(&j) {
                    continue; // Already in a cluster
                }

                let similarity = calculate_jaccard_similarity(
                    &records_with_signatures[i].1,
                    &records_with_signatures[j].1,
                );

                if similarity >= threshold {
                    current_cluster.push(j);
                    assigned.insert(j);
                }
            }

            // Only keep clusters with actual duplicates
            if current_cluster.len() > 1 {
                clusters.push(current_cluster);
            }
        }
    }

    Ok(clusters)
}

fn calculate_jaccard_similarity(minimizers1: &[Minimizer], minimizers2: &[Minimizer]) -> f64 {
    if minimizers1.is_empty() || minimizers2.is_empty() {
        return 0.0;
    }

    // Convert minimizers to hash sets for efficient operations
    let set1: HashSet<u64> = minimizers1.iter().map(|m| m.hash).collect();
    let set2: HashSet<u64> = minimizers2.iter().map(|m| m.hash).collect();

    // Calculate Jaccard index: |intersection| / |union|
    let intersection_size = set1.intersection(&set2).count();
    let union_size = set1.union(&set2).count();

    if union_size == 0 {
        0.0
    } else {
        intersection_size as f64 / union_size as f64
    }
}

fn select_best_quality_read(
    cluster: &[usize],
    records_with_signatures: &[(FastqRecord, Vec<Minimizer>)]
) -> Result<usize> {
    let mut best_idx = cluster[0];
    let mut best_quality = calculate_mean_quality(&records_with_signatures[cluster[0]].0)?;

    for &idx in &cluster[1..] {
        let quality = calculate_mean_quality(&records_with_signatures[idx].0)?;
        if quality > best_quality {
            best_quality = quality;
            best_idx = idx;
        }
    }

    Ok(best_idx)
}

fn calculate_mean_quality(record: &FastqRecord) -> Result<f64> {
    let qualities = &record.quality;
    if qualities.is_empty() {
        return Ok(0.0);
    }

    let sum: u32 = qualities.iter().map(|&q| (q - 33) as u32).sum();
    Ok(sum as f64 / qualities.len() as f64)
}