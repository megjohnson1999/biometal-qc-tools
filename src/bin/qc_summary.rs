//! Biometal QC Summary Tool
//!
//! Multi-sample QC reporting and summary generation

use anyhow::Result;
use biometal_qc_tools::contamination::ContaminationReport;
use biometal_qc_tools::reporting::{QcReporter, SampleQcReport};
use biometal_qc_tools::vlp::VlpReport;
use biometal_qc_tools::QcStats;
use clap::{Arg, Command};
use std::collections::HashSet;
use std::fs;
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

    // Load and process QC result files from input directory
    let mut sample_reports = load_sample_reports(&input_dir)?;

    println!("📈 Generating multi-sample QC summary...");
    println!("📂 Found {} samples to process", sample_reports.len());

    // Evaluate each sample and set pass/fail status
    for sample_report in &mut sample_reports {
        sample_report.overall_pass = reporter.evaluate_sample(sample_report);
    }

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

/// Load and combine QC report files from a directory
fn load_sample_reports(input_dir: &PathBuf) -> Result<Vec<SampleQcReport>> {
    let mut sample_reports = Vec::new();
    let mut sample_names = HashSet::new();

    // First pass: collect all unique sample names
    let entries = fs::read_dir(input_dir)?;
    for entry in entries {
        let entry = entry?;
        let path = entry.path();

        // Skip non-JSON files
        if !path.is_file() || path.extension().map_or(true, |ext| ext != "json") {
            continue;
        }

        let filename = path.file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("unknown");

        if let Some(sample_name) = extract_sample_name(filename) {
            sample_names.insert(sample_name);
        }
    }

    // Second pass: process each unique sample
    for sample_name in sample_names {
        let quality_stats = load_quality_stats(input_dir, &sample_name);
        let contamination_report = load_contamination_report(input_dir, &sample_name);
        let vlp_report = load_vlp_report(input_dir, &sample_name);

        // Only create a sample report if we have at least quality stats
        if let Ok(stats) = quality_stats {
            let contamination = contamination_report.unwrap_or_else(|_| {
                // Create default contamination report if not found
                ContaminationReport {
                    sample_name: sample_name.clone(),
                    total_reads: stats.total_reads,
                    phix_reads: 0,
                    vector_reads: 0,
                    phix_percentage: 0.0,
                    vector_percentage: 0.0,
                }
            });

            let vlp = vlp_report.unwrap_or_else(|_| {
                // Create default VLP report if not found
                VlpReport {
                    sample_name: sample_name.clone(),
                    total_reads: stats.total_reads,
                    gc_distribution_score: 0.0,
                    complexity_diversity: 0.0,
                    compositional_evenness: 0.0,
                    vlp_success_score: 0.0,
                }
            });

            // Create combined sample report
            let sample_report = SampleQcReport {
                quality_stats: stats,
                contamination_report: contamination,
                vlp_report: vlp,
                overall_pass: false, // Will be determined by QcReporter
            };

            sample_reports.push(sample_report);
        }
    }

    Ok(sample_reports)
}

/// Extract sample name from QC report filenames
fn extract_sample_name(filename: &str) -> Option<String> {
    // Remove common QC report suffixes to get sample name
    let name = filename.strip_suffix(".json").unwrap_or(filename);

    if let Some(base) = name.strip_suffix("_quality_stats") {
        return Some(base.to_string());
    }
    if let Some(base) = name.strip_suffix("_contamination_report") {
        return Some(base.to_string());
    }
    if let Some(base) = name.strip_suffix("_vlp_assessment") {
        return Some(base.to_string());
    }
    if let Some(base) = name.strip_suffix("_contamination") {
        return Some(base.to_string());
    }
    if let Some(base) = name.strip_suffix("_vlp") {
        return Some(base.to_string());
    }
    if let Some(base) = name.strip_suffix("_qc") {
        return Some(base.to_string());
    }

    // For files that don't match patterns, use the full name
    Some(name.to_string())
}

/// Load quality statistics for a sample
fn load_quality_stats(input_dir: &PathBuf, sample_name: &str) -> Result<QcStats> {
    let patterns = [
        format!("{}_quality_stats.json", sample_name),
        format!("{}_qc.json", sample_name),
        format!("quality_stats.json"), // Default filename
    ];

    for pattern in &patterns {
        let path = input_dir.join(pattern);
        if path.exists() {
            let content = fs::read_to_string(&path)?;
            return Ok(serde_json::from_str(&content)?);
        }
    }

    Err(anyhow::anyhow!("Quality stats not found for sample: {}", sample_name))
}

/// Load contamination report for a sample
fn load_contamination_report(input_dir: &PathBuf, sample_name: &str) -> Result<ContaminationReport> {
    let patterns = [
        format!("{}_contamination_report.json", sample_name),
        format!("{}_contamination.json", sample_name),
        format!("contamination_report.json"), // Default filename
    ];

    for pattern in &patterns {
        let path = input_dir.join(pattern);
        if path.exists() {
            let content = fs::read_to_string(&path)?;
            return Ok(serde_json::from_str(&content)?);
        }
    }

    Err(anyhow::anyhow!("Contamination report not found for sample: {}", sample_name))
}

/// Load VLP report for a sample
fn load_vlp_report(input_dir: &PathBuf, sample_name: &str) -> Result<VlpReport> {
    let patterns = [
        format!("{}_vlp_assessment.json", sample_name),
        format!("{}_vlp.json", sample_name),
        format!("vlp_assessment.json"), // Default filename
    ];

    for pattern in &patterns {
        let path = input_dir.join(pattern);
        if path.exists() {
            let content = fs::read_to_string(&path)?;
            return Ok(serde_json::from_str(&content)?);
        }
    }

    Err(anyhow::anyhow!("VLP report not found for sample: {}", sample_name))
}