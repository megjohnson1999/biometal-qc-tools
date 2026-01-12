//! Biometal Contamination Screening Tool
//!
//! PhiX and vector detection using biometal pattern matching

use anyhow::Result;
use biometal_qc_tools::contamination::ContaminationScreener;
use clap::{Arg, Command};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-contamination-screen")
        .version("0.1.0")
        .about("Fast contamination screening for FASTQ files")
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
                .help("Output JSON file for contamination report")
                .default_value("contamination_report.json"),
        )
        .arg(
            Arg::new("phix_threshold")
                .long("phix-threshold")
                .value_name("PERCENT")
                .help("PhiX contamination threshold (%)")
                .default_value("0.1"),
        )
        .arg(
            Arg::new("vector_threshold")
                .long("vector-threshold")
                .value_name("PERCENT")
                .help("Vector contamination threshold (%)")
                .default_value("0.05"),
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
    let phix_threshold: f64 = matches
        .get_one::<String>("phix_threshold")
        .unwrap()
        .parse()?;
    let vector_threshold: f64 = matches
        .get_one::<String>("vector_threshold")
        .unwrap()
        .parse()?;
    let min_length: usize = matches
        .get_one::<String>("min_length")
        .unwrap()
        .parse()?;

    println!("🔍 Biometal Contamination Screening Tool");
    println!("Input: {}", input_file.display());
    println!("Output: {}", output_file.display());
    println!("Min Length: {}", min_length);

    // Validate input file exists
    if !input_file.exists() {
        anyhow::bail!("Input file does not exist: {}", input_file.display());
    }

    // Create contamination screener
    let screener = ContaminationScreener::new(phix_threshold, vector_threshold, min_length);

    // Screen for contamination
    println!("🦠 Screening for contamination...");
    let report = screener.screen_fastq(&input_file)?;

    // Check if contamination is acceptable
    let acceptable = screener.is_contamination_acceptable(&report);

    // Output results
    let json_output = serde_json::to_string_pretty(&report)?;
    std::fs::write(&output_file, &json_output)?;

    println!("✅ Contamination screening complete!");
    println!("📈 Sample: {}", report.sample_name);
    println!("🦠 PhiX contamination: {:.3}%", report.phix_percentage);
    println!("🧬 Vector contamination: {:.3}%", report.vector_percentage);
    println!(
        "{}",
        if acceptable {
            "✅ Contamination levels acceptable"
        } else {
            "⚠️ Contamination levels exceed thresholds"
        }
    );
    println!("💾 Results saved to: {}", output_file.display());

    Ok(())
}