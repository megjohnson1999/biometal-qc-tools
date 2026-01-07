//! Biometal Quality Stats Tool
//!
//! Fast quality assessment tool to replace FastQC functionality
//! Uses proven biometal primitives: base_counting, gc_content, quality_filter, complexity

use anyhow::Result;
use biometal_qc_tools::quality::QualityAnalyzer;
use clap::{Arg, Command};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-quality-stats")
        .version("0.1.0")
        .about("Fast quality assessment for FASTQ files using biometal primitives")
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
                .value_name("JSON")
                .help("Output JSON file for statistics")
                .default_value("quality_stats.json"),
        )
        .arg(
            Arg::new("min_quality")
                .short('q')
                .long("min-quality")
                .value_name("QUALITY")
                .help("Minimum quality score threshold")
                .default_value("20"),
        )
        .arg(
            Arg::new("min_length")
                .short('l')
                .long("min-length")
                .value_name("LENGTH")
                .help("Minimum read length")
                .default_value("50"),
        )
        .get_matches();

    // Parse arguments
    let input_file = PathBuf::from(matches.get_one::<String>("input").unwrap());
    let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
    let min_quality: u8 = matches
        .get_one::<String>("min_quality")
        .unwrap()
        .parse()?;
    let min_length: usize = matches
        .get_one::<String>("min_length")
        .unwrap()
        .parse()?;

    println!("🧬 Biometal Quality Stats Tool");
    println!("Input: {}", input_file.display());
    println!("Output: {}", output_file.display());
    println!("Min Quality: {}, Min Length: {}", min_quality, min_length);

    // Validate input file exists
    if !input_file.exists() {
        anyhow::bail!("Input file does not exist: {}", input_file.display());
    }

    // Create quality analyzer
    let analyzer = QualityAnalyzer::new(min_quality, min_length);

    // Analyze the FASTQ file
    println!("📊 Analyzing quality statistics...");
    let stats = analyzer.analyze_fastq(&input_file)?;

    // Output results
    let json_output = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&output_file, &json_output)?;

    println!("✅ Quality analysis complete!");
    println!("📈 Sample: {}", stats.sample_name);
    println!("📚 Total reads: {}", stats.total_reads);
    println!("🧬 Total bases: {}", stats.total_bases);
    println!("🔬 GC content: {:.2}%", stats.gc_content);
    println!("⭐ Mean quality: {:.2}", stats.mean_quality);
    println!("🎯 Q30 bases: {:.2}%", stats.q30_bases);
    println!("🌀 Complexity: {:.2}", stats.complexity_score);
    println!("💾 Results saved to: {}", output_file.display());

    Ok(())
}