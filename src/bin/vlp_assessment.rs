//! Biometal VLP Assessment Tool
//!
//! Composition-based VLP success metrics using biometal primitives

use anyhow::Result;
use biometal_qc_tools::vlp::VlpAssessor;
use clap::{Arg, Command};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-vlp-assessment")
        .version("0.1.0")
        .about("VLP success assessment using composition-based metrics")
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
                .help("Output JSON file for VLP assessment")
                .default_value("vlp_assessment.json"),
        )
        .arg(
            Arg::new("min_complexity")
                .long("min-complexity")
                .value_name("SCORE")
                .help("Minimum complexity score threshold")
                .default_value("0.7"),
        )
        .arg(
            Arg::new("gc_min")
                .long("gc-min")
                .value_name("PERCENT")
                .help("Minimum GC content for optimal range")
                .default_value("0.35"),
        )
        .arg(
            Arg::new("gc_max")
                .long("gc-max")
                .value_name("PERCENT")
                .help("Maximum GC content for optimal range")
                .default_value("0.65"),
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
    let min_complexity: f64 = matches
        .get_one::<String>("min_complexity")
        .unwrap()
        .parse()?;
    let gc_min: f64 = matches.get_one::<String>("gc_min").unwrap().parse()?;
    let gc_max: f64 = matches.get_one::<String>("gc_max").unwrap().parse()?;
    let min_length: usize = matches
        .get_one::<String>("min_length")
        .unwrap()
        .parse()?;

    println!("🦠 Biometal VLP Assessment Tool");
    println!("Input: {}", input_file.display());
    println!("Output: {}", output_file.display());
    println!("Min Length: {}", min_length);

    // Validate input file exists
    if !input_file.exists() {
        anyhow::bail!("Input file does not exist: {}", input_file.display());
    }

    // Create VLP assessor
    let assessor = VlpAssessor::new(min_complexity, (gc_min, gc_max), min_length);

    // Assess VLP success
    println!("🧬 Assessing VLP success metrics...");
    let report = assessor.assess_vlp(&input_file)?;

    // Check if VLP preparation was successful
    let successful = assessor.is_vlp_successful(&report);

    // Output results
    let json_output = serde_json::to_string_pretty(&report)?;
    std::fs::write(&output_file, &json_output)?;

    println!("✅ VLP assessment complete!");
    println!("📈 Sample: {}", report.sample_name);
    println!("🔬 GC distribution score: {:.3}", report.gc_distribution_score);
    println!("🌀 Complexity diversity: {:.3}", report.complexity_diversity);
    println!("⚖️ Compositional evenness: {:.3}", report.compositional_evenness);
    println!("🎯 Overall VLP score: {:.3}", report.vlp_success_score);
    println!(
        "{}",
        if successful {
            "✅ VLP preparation successful"
        } else {
            "⚠️ VLP preparation may need optimization"
        }
    );
    println!("💾 Results saved to: {}", output_file.display());

    Ok(())
}