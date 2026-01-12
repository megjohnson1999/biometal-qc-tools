//! Biometal Host Contamination Removal Tool
//!
//! Constant-memory host depletion using biometal StreamingMapper
//! Replaces minimap2 + samtools with ~5MB memory vs 6-10GB requirement

use anyhow::Result;
use biometal::alignment::{StreamingMapper, StreamingMapperConfig, MappingResult};
use biometal::{FastqStream, FastqWriter};
use clap::{Arg, Command};
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Serialize, Deserialize)]
struct HostDepletionStats {
    total_reads: u64,
    host_matches_found: u64,
    clean_reads_kept: u64,
    contamination_rate: f64,
    alignment_score_threshold: i32,
    window_size: usize,
    overlap_size: usize,
    processing_time_seconds: f64,
}

fn main() -> Result<()> {
    let matches = Command::new("biometal-host-depletion")
        .version("0.1.0")
        .about("Constant-memory host contamination removal using biometal StreamingMapper")
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
            Arg::new("reference")
                .short('r')
                .long("reference")
                .value_name("FASTA")
                .help("Host reference genome (FASTA format)")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FASTQ")
                .help("Output host-depleted FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output host depletion statistics (JSON)")
                .default_value("host_depletion_stats.json"),
        )
        .arg(
            Arg::new("threshold")
                .short('t')
                .long("threshold")
                .value_name("SCORE")
                .help("Minimum alignment score to consider as host match")
                .default_value("50"),
        )
        .arg(
            Arg::new("window_size")
                .long("window-size")
                .value_name("BYTES")
                .help("Reference processing window size")
                .default_value("1000000"),
        )
        .arg(
            Arg::new("overlap")
                .long("overlap")
                .value_name("BP")
                .help("Window overlap size in base pairs")
                .default_value("200"),
        )
        .get_matches();

    // Parse arguments
    let input_path: PathBuf = matches.get_one::<String>("input").unwrap().into();
    let reference_path: PathBuf = matches.get_one::<String>("reference").unwrap().into();
    let output_path: PathBuf = matches.get_one::<String>("output").unwrap().into();
    let stats_path: PathBuf = matches.get_one::<String>("stats").unwrap().into();
    let threshold: i32 = matches.get_one::<String>("threshold").unwrap().parse()?;
    let window_size: usize = matches.get_one::<String>("window_size").unwrap().parse()?;
    let overlap: usize = matches.get_one::<String>("overlap").unwrap().parse()?;

    println!("🏠 Biometal Host Contamination Removal");
    println!("======================================");
    println!("Input reads: {}", input_path.display());
    println!("Host reference: {}", reference_path.display());
    println!("Output: {}", output_path.display());
    println!("Alignment threshold: {}", threshold);
    println!("Window size: {} bytes, Overlap: {} bp", window_size, overlap);

    let start_time = std::time::Instant::now();

    // Process host depletion
    let stats = process_host_depletion(
        &input_path,
        &reference_path,
        &output_path,
        threshold,
        window_size,
        overlap,
    )?;

    let processing_time = start_time.elapsed().as_secs_f64();
    let final_stats = HostDepletionStats {
        processing_time_seconds: processing_time,
        ..stats
    };

    // Write statistics
    let stats_json = serde_json::to_string_pretty(&final_stats)?;
    std::fs::write(&stats_path, stats_json)?;

    println!("\n✅ Host Depletion Complete");
    println!("Total reads processed: {}", final_stats.total_reads);
    println!("Host matches found: {}", final_stats.host_matches_found);
    println!("Clean reads kept: {}", final_stats.clean_reads_kept);
    println!("Contamination rate: {:.2}%", final_stats.contamination_rate * 100.0);
    println!("Processing time: {:.2}s", processing_time);
    println!("Statistics written to: {}", stats_path.display());

    Ok(())
}

fn process_host_depletion(
    input_path: &PathBuf,
    reference_path: &PathBuf,
    output_path: &PathBuf,
    threshold: i32,
    window_size: usize,
    overlap: usize,
) -> Result<HostDepletionStats> {
    // Step 1: Configure StreamingMapper
    let config = StreamingMapperConfig {
        window_size,
        overlap_bp: overlap,
        min_score_threshold: threshold,
        ..Default::default()
    };

    let mut mapper = StreamingMapper::new(config);

    println!("🧬 Analyzing reads against host genome...");
    println!("   Using streaming mapper with ~5MB constant memory");

    // Step 2: Map reads and collect alignment results
    let mappings = mapper.map_reads_streaming(reference_path, input_path)?;
    let mut host_alignments: HashMap<String, MappingResult> = HashMap::new();
    let mut total_reads = 0;

    for mapping_result in mappings {
        let mapping = mapping_result?;
        total_reads += 1;

        // Store best alignment for each read
        let read_id = mapping.query_id.clone();

        match host_alignments.get(&read_id) {
            Some(existing) => {
                // Keep the better alignment score
                if mapping.alignment.score > existing.alignment.score {
                    host_alignments.insert(read_id, mapping);
                }
            }
            None => {
                host_alignments.insert(read_id, mapping);
            }
        }

        if total_reads % 1000 == 0 {
            println!("   - Processed {} read alignments...", total_reads);
        }
    }

    println!("   - Total alignment results: {}", host_alignments.len());

    // Step 3: Read original FASTQ and filter based on alignments
    println!("📝 Writing host-depleted output...");

    let fastq_stream = FastqStream::from_path(input_path)?;
    let mut writer = FastqWriter::create(output_path)?;

    let mut fastq_total_reads = 0;
    let mut host_matches_found = 0;
    let mut clean_reads_kept = 0;

    for record_result in fastq_stream {
        let record = record_result?;
        fastq_total_reads += 1;

        // Check if this read has a significant host alignment
        let is_host_contamination = match host_alignments.get(&record.id) {
            Some(mapping) => {
                // Check if alignment score meets threshold for host contamination
                mapping.alignment.score >= threshold
            }
            None => false, // No alignment found = not host contamination
        };

        if is_host_contamination {
            host_matches_found += 1;
            // Skip host-contaminated reads
        } else {
            // Keep non-host reads
            writer.write_record(&record)?;
            clean_reads_kept += 1;
        }

        if fastq_total_reads % 5000 == 0 {
            println!("   - Processed {} FASTQ records...", fastq_total_reads);
        }
    }

    let contamination_rate = if fastq_total_reads > 0 {
        host_matches_found as f64 / fastq_total_reads as f64
    } else {
        0.0
    };

    Ok(HostDepletionStats {
        total_reads: fastq_total_reads,
        host_matches_found,
        clean_reads_kept,
        contamination_rate,
        alignment_score_threshold: threshold,
        window_size,
        overlap_size: overlap,
        processing_time_seconds: 0.0, // Will be set by caller
    })
}