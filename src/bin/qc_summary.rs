//! Biometal QC Summary Tool
//!
//! Multi-sample QC reporting and summary generation

use anyhow::Result;
use biometal_qc_tools::reporting::{QcReporter, SampleQcReport};
use clap::{Arg, Command};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-qc-summary")
        .version("0.1.0")
        .about("Multi-sample QC summary and reporting")
        .author("Megan Johnson")
        .arg(
            Arg::new("input_dir")
                .short('i')
                .long("input-dir")
                .value_name("DIRECTORY")
                .help("Directory containing QC result JSON files")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("JSON")
                .help("Output JSON file for summary report")
                .default_value("qc_summary.json"),
        )
        .arg(
            Arg::new("quality_threshold")
                .long("quality-threshold")
                .value_name("SCORE")
                .help("Quality score threshold for pass/fail")
                .default_value("25.0"),
        )
        .arg(
            Arg::new("contamination_threshold")
                .long("contamination-threshold")
                .value_name("PERCENT")
                .help("Contamination threshold for pass/fail")
                .default_value("0.1"),
        )
        .get_matches();

    // Parse arguments
    let input_dir = PathBuf::from(matches.get_one::<String>("input_dir").unwrap());
    let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
    let quality_threshold: f64 = matches
        .get_one::<String>("quality_threshold")
        .unwrap()
        .parse()?;
    let contamination_threshold: f64 = matches
        .get_one::<String>("contamination_threshold")
        .unwrap()
        .parse()?;

    println!("📊 Biometal QC Summary Tool");
    println!("Input directory: {}", input_dir.display());
    println!("Output: {}", output_file.display());

    // Validate input directory exists
    if !input_dir.exists() || !input_dir.is_dir() {
        anyhow::bail!("Input directory does not exist: {}", input_dir.display());
    }

    // Create QC reporter
    let reporter = QcReporter::new(quality_threshold, contamination_threshold);

    // TODO: Load and process QC result files from input directory
    // For now, create an empty sample list
    let sample_reports: Vec<SampleQcReport> = Vec::new();

    println!("📈 Generating multi-sample QC summary...");

    // Generate comprehensive report
    let multi_sample_report = reporter.generate_report(sample_reports);

    // Export to JSON
    reporter.export_json(&multi_sample_report, &output_file)?;

    println!("✅ QC summary complete!");
    println!("📊 Summary Statistics:");
    println!("  Total samples: {}", multi_sample_report.summary.total_samples);
    println!("  Passed samples: {}", multi_sample_report.summary.passed_samples);
    println!("  Failed samples: {}", multi_sample_report.summary.failed_samples);
    println!("  Pass rate: {:.1}%", multi_sample_report.summary.pass_rate);
    println!("  Average quality: {:.2}", multi_sample_report.summary.average_quality);
    println!("  Average GC: {:.2}%", multi_sample_report.summary.average_gc_content);
    println!("💾 Summary saved to: {}", output_file.display());

    Ok(())
}