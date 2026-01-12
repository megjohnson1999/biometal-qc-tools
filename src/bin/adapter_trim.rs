//! Biometal Adapter Trimming Tool
//!
//! Fast adapter detection and trimming tool to replace fastp's adapter trimming functionality
//! Uses proven biometal primitives: AdapterDetector, pattern matching, and trimming operations

use anyhow::Result;
use biometal_qc_tools::adapters::AdapterTrimmer;
use clap::{Arg, Command};
use serde_json;
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-adapter-trim")
        .version("0.1.0")
        .about("Fast adapter trimming for FASTQ files using biometal primitives")
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
                .help("Output trimmed FASTQ file")
                .required(false),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output adapter trimming statistics (JSON)")
                .default_value("adapter_stats.json"),
        )
        .arg(
            Arg::new("min_adapter_length")
                .long("min-adapter-length")
                .value_name("LENGTH")
                .help("Minimum adapter match length to consider for trimming")
                .default_value("8"),
        )
        .arg(
            Arg::new("min_overlap")
                .long("min-overlap")
                .value_name("OVERLAP")
                .help("Minimum overlap length to trigger trimming")
                .default_value("5"),
        )
        .arg(
            Arg::new("trim_3_only")
                .long("trim-3-only")
                .help("Only trim 3' end adapters (default: trim both ends)")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .help("Verbose output")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    // Parse arguments
    let input_path = PathBuf::from(matches.get_one::<String>("input").unwrap());
    let output_path = matches.get_one::<String>("output").map(PathBuf::from);
    let stats_path = PathBuf::from(matches.get_one::<String>("stats").unwrap());
    let min_adapter_length: usize = matches.get_one::<String>("min_adapter_length").unwrap().parse()?;
    let min_overlap: usize = matches.get_one::<String>("min_overlap").unwrap().parse()?;
    let trim_both_ends = !matches.get_flag("trim_3_only");
    let verbose = matches.get_flag("verbose");

    if verbose {
        println!("Biometal Adapter Trimming Tool v0.1.0");
        println!("=====================================");
        println!("Input file: {}", input_path.display());
        if let Some(ref out_path) = output_path {
            println!("Output file: {}", out_path.display());
        } else {
            println!("Output file: None (stats only)");
        }
        println!("Stats file: {}", stats_path.display());
        println!("Min adapter length: {}", min_adapter_length);
        println!("Min overlap: {}", min_overlap);
        println!("Trim both ends: {}", trim_both_ends);
        println!();
    }

    // Validate input file exists
    if !input_path.exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {}", input_path.display()));
    }

    // Create adapter trimmer
    let trimmer = AdapterTrimmer::new(min_adapter_length, min_overlap, trim_both_ends);

    if verbose {
        println!("Processing FASTQ file...");
    }

    // Process the FASTQ file
    let stats = trimmer.process_fastq(&input_path, output_path.as_ref())?;

    if verbose {
        println!("Adapter trimming completed!");
        println!();
        println!("Summary Statistics:");
        println!("==================");
        println!("Total reads processed: {}", stats.total_reads);
        println!("Reads with adapters: {} ({:.1}%)",
                 stats.reads_with_adapters,
                 100.0 * stats.reads_with_adapters as f64 / stats.total_reads as f64);
        println!("Total bases trimmed: {}", stats.total_bases_trimmed);

        if stats.reads_with_adapters > 0 {
            println!("Average bases trimmed per affected read: {:.1}", stats.average_trim_length);
        }

        if !stats.adapters_found.is_empty() {
            println!();
            println!("Adapter Types Found:");
            for (adapter_name, count) in &stats.adapters_found {
                println!("  {}: {} occurrences", adapter_name, count);
            }
        }
        println!();
    }

    // Write statistics to JSON file
    let stats_json = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_path, stats_json)?;

    if verbose {
        println!("Statistics written to: {}", stats_path.display());
    }

    // Summary message
    if stats.reads_with_adapters > 0 {
        println!("✅ Adapter trimming completed: {} reads processed, {} adapters trimmed",
                 stats.total_reads, stats.reads_with_adapters);
    } else {
        println!("✅ No adapters found in {} reads", stats.total_reads);
    }

    Ok(())
}