//! Biometal rRNA Removal Tool
//!
//! Advanced rRNA detection and removal using biometal's algorithmic primitives
//! Demonstrates superior approach over traditional tools like BBDuk:
//!
//! 1. Minimizer-based database indexing for fast screening
//! 2. Smith-Waterman alignment for sensitive confirmation
//! 3. K-mer spectrum analysis for content assessment
//! 4. NEON-optimized operations with memory-efficient streaming

use anyhow::Result;
use biometal_qc_tools::rrna::RrnaRemover;
use clap::{Arg, Command};
use serde_json;
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-rrna-remove")
        .version("0.1.0")
        .about("Advanced rRNA detection and removal using biometal algorithmic primitives")
        .long_about("Showcases biometal's algorithmic advantages over traditional tools:\\n\\\n                     • Minimizer-based rRNA database fingerprinting for fast screening\\n\\\n                     • Smith-Waterman alignment for sensitive rRNA detection with mismatches\\n\\\n                     • K-mer spectrum analysis for rRNA content assessment\\n\\\n                     • NEON-optimized operations with streaming database processing\\n\\\n                     \\n\\\n                     Unlike BBDuk's rigid k-mer matching, provides superior sensitivity\\n\\\n                     and memory efficiency for massive Silva databases.")
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
            Arg::new("database")
                .short('d')
                .long("database")
                .value_name("FASTA")
                .help("rRNA reference database (FASTA format, e.g., Silva SSU/LSU)")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FASTQ")
                .help("Output rRNA-filtered FASTQ file")
                .required(false),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output rRNA removal statistics (JSON)")
                .default_value("rrna_removal_stats.json"),
        )
        .arg(
            Arg::new("minimizer_length")
                .long("minimizer-length")
                .value_name("LENGTH")
                .help("Minimizer length for fast screening (default: 15)")
                .default_value("15"),
        )
        .arg(
            Arg::new("alignment_threshold")
                .long("alignment-threshold")
                .value_name("SCORE")
                .help("Smith-Waterman alignment score threshold (0.0-1.0, default: 0.8)")
                .default_value("0.8"),
        )
        .arg(
            Arg::new("kmer_size")
                .long("kmer-size")
                .value_name("SIZE")
                .help("K-mer size for content analysis (default: 21)")
                .default_value("21"),
        )
        .arg(
            Arg::new("rrna_threshold")
                .long("rrna-threshold")
                .value_name("PERCENT")
                .help("rRNA content threshold for flagging samples (default: 10.0%)")
                .default_value("10.0"),
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
    let database_path = PathBuf::from(matches.get_one::<String>("database").unwrap());
    let output_path = matches.get_one::<String>("output").map(PathBuf::from);
    let stats_path = PathBuf::from(matches.get_one::<String>("stats").unwrap());
    let minimizer_length: usize = matches.get_one::<String>("minimizer_length").unwrap().parse()?;
    let alignment_threshold: f64 = matches.get_one::<String>("alignment_threshold").unwrap().parse()?;
    let kmer_size: usize = matches.get_one::<String>("kmer_size").unwrap().parse()?;
    let rrna_threshold: f64 = matches.get_one::<String>("rrna_threshold").unwrap().parse()?;
    let verbose = matches.get_flag("verbose");

    if verbose {
        println!("🧬 Biometal rRNA Removal Tool v0.1.0");
        println!("=====================================");
        println!("Algorithmic Approach:");
        println!("  1. Minimizer indexing for fast screening");
        println!("  2. Smith-Waterman alignment for confirmation");
        println!("  3. K-mer spectrum analysis for content assessment");
        println!("  4. NEON-optimized + streaming for efficiency");
        println!();
        println!("Configuration:");
        println!("  Input file: {}", input_path.display());
        println!("  rRNA database: {}", database_path.display());
        if let Some(ref out_path) = output_path {
            println!("  Output file: {}", out_path.display());
        } else {
            println!("  Output file: None (statistics only)");
        }
        println!("  Stats file: {}", stats_path.display());
        println!("  Minimizer length: {}", minimizer_length);
        println!("  Alignment threshold: {:.2}", alignment_threshold);
        println!("  K-mer size: {}", kmer_size);
        println!("  rRNA content threshold: {:.1}%", rrna_threshold);
        println!();
    }

    // Validate input files
    if !input_path.exists() {
        return Err(anyhow::anyhow!("Input FASTQ file does not exist: {}", input_path.display()));
    }

    if !database_path.exists() {
        return Err(anyhow::anyhow!("rRNA database file does not exist: {}", database_path.display()));
    }

    // Validate parameters
    if alignment_threshold < 0.0 || alignment_threshold > 1.0 {
        return Err(anyhow::anyhow!("Alignment threshold must be between 0.0 and 1.0, got: {}", alignment_threshold));
    }

    if minimizer_length < 10 || minimizer_length > 25 {
        if verbose {
            println!("⚠️  Warning: Unusual minimizer length ({}). Recommended range: 10-25", minimizer_length);
        }
    }

    if kmer_size < 15 || kmer_size > 31 {
        if verbose {
            println!("⚠️  Warning: Unusual k-mer size ({}). Recommended range: 15-31", kmer_size);
        }
    }

    // Create rRNA remover with biometal algorithms
    let mut remover = RrnaRemover::new(minimizer_length, alignment_threshold, kmer_size);
    remover.rrna_content_threshold = rrna_threshold / 100.0; // Convert percentage to fraction

    if verbose {
        println!("🚀 Starting biometal rRNA removal pipeline...");
        println!("   This showcases biometal's algorithmic advantages:");
        println!("   • {}× faster minimizer extraction (NEON-optimized)", if cfg!(target_arch = "aarch64") { "8-15" } else { "2-4" });
        println!("   • Sensitive Smith-Waterman alignment vs exact matching");
        println!("   • Memory-efficient streaming through large databases");
        println!("   • Advanced k-mer spectrum analysis");
    }

    // Process the FASTQ file with advanced biometal algorithms
    let stats = remover.process_fastq(&input_path, &database_path, output_path.as_ref())?;

    if verbose {
        println!("✅ rRNA removal pipeline completed!");
        println!();
        println!("📊 Algorithmic Performance Summary:");
        println!("===================================");
        println!("Database sequences indexed: {}", stats.database_sequences_processed);
        println!("Total reads processed: {}", stats.total_reads);
        println!("Minimizer screening hits: {}", stats.minimizer_matches);
        println!("Smith-Waterman confirmations: {}", stats.alignment_confirmations);
        println!("rRNA reads detected: {} ({:.1}%)", stats.rrna_reads_detected, stats.rrna_detection_rate);
        println!("rRNA reads removed: {}", stats.rrna_reads_removed);
        println!("Clean reads retained: {}", stats.total_reads - stats.rrna_reads_removed);
        println!();
        println!("📈 Content Analysis:");
        println!("===================");
        println!("Sample rRNA content score: {:.3}", stats.rrna_content_score);
        println!("rRNA screening efficiency: {:.1}% (minimizer hits / total reads)",
                 100.0 * stats.minimizer_matches as f64 / stats.total_reads as f64);
        println!("Alignment confirmation rate: {:.1}% (confirmations / minimizer hits)",
                 if stats.minimizer_matches > 0 {
                     100.0 * stats.alignment_confirmations as f64 / stats.minimizer_matches as f64
                 } else {
                     0.0
                 });

        println!();

        // Sample quality assessment
        if remover.is_rrna_content_high(&stats) {
            println!("⚠️  High rRNA content detected ({:.1}% > {:.1}%)",
                     stats.rrna_content_score * 100.0, rrna_threshold);
            println!("   This may indicate:");
            println!("   • Incomplete rRNA depletion during library prep");
            println!("   • Environmental/microbial contamination");
            println!("   • Non-poly(A) selected RNA library");
        } else {
            println!("✅ rRNA content acceptable ({:.1}% ≤ {:.1}%)",
                     stats.rrna_content_score * 100.0, rrna_threshold);
        }

        println!();
        println!("🔬 Biometal Advantages Demonstrated:");
        println!("===================================");
        println!("• Minimizer indexing: Fast database screening vs linear search");
        println!("• Smith-Waterman: Sensitive alignment vs rigid k-mer matching");
        println!("• Streaming I/O: Constant memory usage for massive Silva databases");
        println!("• NEON optimization: Hardware acceleration on ARM platforms");
        println!("• K-mer analysis: Advanced content assessment vs simple counting");
    }

    // Write statistics to JSON file
    let stats_json = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_path, stats_json)?;

    if verbose {
        println!("💾 Detailed statistics saved to: {}", stats_path.display());
    }

    // Summary message
    if stats.rrna_reads_detected > 0 {
        println!("🧬 Biometal rRNA removal completed: {} reads processed, {} rRNA sequences removed ({:.1}%)",
                 stats.total_reads, stats.rrna_reads_removed, stats.rrna_detection_rate);
    } else {
        println!("🧬 No rRNA sequences detected in {} reads (clean sample)", stats.total_reads);
    }

    // Quality warning
    if remover.is_rrna_content_high(&stats) {
        println!("⚠️  Warning: High rRNA content detected ({:.1}%) - review library preparation",
                 stats.rrna_content_score * 100.0);
    }

    Ok(())
}