//! Biometal Primer Removal Tool
//!
//! Fast primer B removal for FASTQ files using biometal primitives
//! Replicates BBDuk's two-step primer removal process from lab-virome-QC pipeline

use anyhow::Result;
use biometal_qc_tools::primers::PrimerRemover;
use clap::{Arg, Command};
use serde_json;
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-primer-remove")
        .version("0.1.0")
        .about("Fast primer B removal for FASTQ files using biometal primitives")
        .long_about("Replicates BBDuk's two-step primer removal process:\n\
                     1. Remove forward primer B sequences from 5' end (ktrim=l)\n\
                     2. Remove reverse complement primer B sequences from 3' end (ktrim=r)\n\
                     \n\
                     Uses the same 24 Primer B variants and k-mer matching (k=16, mink=9) as the lab-virome-QC pipeline.")
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
                .help("Output primer-trimmed FASTQ file")
                .required(false),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output primer removal statistics (JSON)")
                .default_value("primer_removal_stats.json"),
        )
        .arg(
            Arg::new("min_match_length")
                .long("min-match-length")
                .value_name("LENGTH")
                .help("Minimum k-mer match length (BBDuk mink parameter)")
                .default_value("9"),
        )
        .arg(
            Arg::new("max_match_length")
                .long("max-match-length")
                .value_name("LENGTH")
                .help("Maximum k-mer match length (BBDuk k parameter)")
                .default_value("16"),
        )
        .arg(
            Arg::new("contamination_threshold")
                .long("contamination-threshold")
                .value_name("PERCENT")
                .help("Cross-contamination threshold for flagging (%)")
                .default_value("5.0"),
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
    let min_match_length: usize = matches.get_one::<String>("min_match_length").unwrap().parse()?;
    let max_match_length: usize = matches.get_one::<String>("max_match_length").unwrap().parse()?;
    let contamination_threshold: f64 = matches.get_one::<String>("contamination_threshold").unwrap().parse()?;
    let verbose = matches.get_flag("verbose");

    if verbose {
        println!("🧬 Biometal Primer B Removal Tool v0.1.0");
        println!("==========================================");
        println!("Input file: {}", input_path.display());
        if let Some(ref out_path) = output_path {
            println!("Output file: {}", out_path.display());
        } else {
            println!("Output file: None (stats only)");
        }
        println!("Stats file: {}", stats_path.display());
        println!("K-mer range: {} to {} (BBDuk: mink={}, k={})", min_match_length, max_match_length, min_match_length, max_match_length);
        println!("Contamination threshold: {:.1}%", contamination_threshold);
        println!();
    }

    // Validate input file exists
    if !input_path.exists() {
        return Err(anyhow::anyhow!("Input file does not exist: {}", input_path.display()));
    }

    // Validate parameters
    if min_match_length > max_match_length {
        return Err(anyhow::anyhow!("min-match-length ({}) cannot be greater than max-match-length ({})", min_match_length, max_match_length));
    }

    if max_match_length != 16 || min_match_length != 9 {
        if verbose {
            println!("⚠️  Warning: Using non-standard k-mer parameters. BBDuk pipeline uses k=16, mink=9");
        }
    }

    // Create primer remover
    let remover = PrimerRemover::new(min_match_length, max_match_length, contamination_threshold);

    if verbose {
        println!("🔍 Processing FASTQ file...");
        println!("   Step 1: Removing forward primers from 5' end");
        println!("   Step 2: Removing reverse complement primers from 3' end");
    }

    // Process the FASTQ file
    let stats = remover.process_fastq(&input_path, output_path.as_ref())?;

    if verbose {
        println!("✅ Primer removal completed!");
        println!();
        println!("📊 Summary Statistics:");
        println!("======================");
        println!("Total reads processed: {}", stats.total_reads);
        println!("Reads with forward primers: {} ({:.1}%)",
                 stats.reads_with_forward_primers,
                 100.0 * stats.reads_with_forward_primers as f64 / stats.total_reads as f64);
        println!("Reads with RC primers: {} ({:.1}%)",
                 stats.reads_with_rc_primers,
                 100.0 * stats.reads_with_rc_primers as f64 / stats.total_reads as f64);
        println!("Total bases trimmed: {}", stats.total_bases_trimmed);
        println!("Cross-contamination level: {:.2}%", stats.contamination_level);

        if !stats.forward_primers_found.is_empty() {
            println!();
            println!("🧬 Forward Primers Found:");
            for (primer_id, count) in &stats.forward_primers_found {
                println!("   {}: {} occurrences", primer_id, count);
            }
        }

        if !stats.rc_primers_found.is_empty() {
            println!();
            println!("🔄 Reverse Complement Primers Found:");
            for (primer_id, count) in &stats.rc_primers_found {
                println!("   {}: {} occurrences", primer_id, count);
            }
        }

        println!();

        // Contamination assessment
        if remover.is_contamination_acceptable(&stats) {
            println!("✅ Cross-contamination levels acceptable (<{:.1}%)", contamination_threshold);
        } else {
            println!("⚠️  High cross-contamination detected ({:.2}% > {:.1}%)", stats.contamination_level, contamination_threshold);
            println!("   This may indicate sample cross-contamination during library prep");
        }
    }

    // Write statistics to JSON file
    let stats_json = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_path, stats_json)?;

    if verbose {
        println!("💾 Statistics saved to: {}", stats_path.display());
    }

    // Summary message
    let total_primer_reads = stats.reads_with_forward_primers + stats.reads_with_rc_primers;
    if total_primer_reads > 0 {
        println!("🧬 Primer removal completed: {} reads processed, {} primer sequences removed",
                 stats.total_reads, total_primer_reads);
    } else {
        println!("🧬 No primer sequences found in {} reads", stats.total_reads);
    }

    // Contamination warning
    if !remover.is_contamination_acceptable(&stats) {
        println!("⚠️  Warning: High cross-contamination detected ({:.2}%)", stats.contamination_level);
    }

    Ok(())
}