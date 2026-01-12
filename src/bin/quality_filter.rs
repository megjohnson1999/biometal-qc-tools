//! Biometal Quality Filter Tool
//!
//! Filter FASTQ reads based on mean quality scores using biometal primitives.

use anyhow::Result;
use biometal::io::{DataSource, FastqStream};
use biometal::operations::mean_quality;
use biometal_qc_tools::{get_file_info, QualityFilterStats};
use clap::{Arg, Command};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-quality-filter")
        .version("0.1.0")
        .about("Filter FASTQ reads based on mean quality scores")
        .author("Megan Johnson")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FASTQ")
                .help("Input FASTQ file (gzip supported)")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FASTQ")
                .help("Output filtered FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("min_quality")
                .long("min-quality")
                .value_name("SCORE")
                .help("Minimum mean quality score threshold")
                .default_value("20.0"),
        )
        .arg(
            Arg::new("stats")
                .long("stats")
                .value_name("JSON")
                .help("Output statistics JSON file")
                .default_value("quality_filter_stats.json"),
        )
        .get_matches();

    // Parse arguments
    let input_file = PathBuf::from(matches.get_one::<String>("input").unwrap());
    let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
    let min_quality: f64 = matches
        .get_one::<String>("min_quality")
        .unwrap()
        .parse()?;
    let stats_file = PathBuf::from(matches.get_one::<String>("stats").unwrap());

    println!("🎯 Biometal Quality Filter");
    println!("Input: {}", input_file.display());
    println!("Output: {}", output_file.display());
    println!("Min quality threshold: {:.1}", min_quality);

    // Get sample name from input file
    let (sample_name, _) = get_file_info(&input_file)?;

    // Create filter and process
    let filter = QualityFilter::new(min_quality);
    let stats = filter.filter_reads(&input_file, &output_file, &sample_name)?;

    // Output statistics
    println!("📊 Quality Filtering Results:");
    println!("  Total reads: {}", stats.total_reads);
    println!("  Reads passed: {} ({:.1}%)",
             stats.reads_passed,
             stats.pass_rate);
    println!("  Reads failed: {} ({:.1}%)",
             stats.reads_failed,
             100.0 - stats.pass_rate);
    println!("  Quality threshold: {:.1}", stats.quality_threshold);

    // Export statistics to JSON
    let json_content = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_file, json_content)?;
    println!("💾 Statistics saved to: {}", stats_file.display());

    Ok(())
}

/// Quality filter implementation
pub struct QualityFilter {
    pub min_quality: f64,
}

impl QualityFilter {
    pub fn new(min_quality: f64) -> Self {
        Self { min_quality }
    }

    /// Filter reads based on mean quality score
    pub fn filter_reads(
        &self,
        input_path: &PathBuf,
        output_path: &PathBuf,
        sample_name: &str,
    ) -> Result<QualityFilterStats> {
        let mut total_reads = 0u64;
        let mut reads_passed = 0u64;

        // Open input stream
        let data_source = DataSource::from_path(input_path);
        let fastq_stream = FastqStream::new(data_source)?;

        // Open output writer
        let output_file = File::create(output_path)?;
        let mut writer = BufWriter::new(output_file);

        // Process each read
        for record_result in fastq_stream {
            let record = record_result?;
            if record.is_empty() {
                continue;
            }

            total_reads += 1;

            // Calculate mean quality using biometal primitive
            let read_mean_quality = mean_quality(&record.quality);

            // Check if read passes quality threshold
            if read_mean_quality >= self.min_quality {
                reads_passed += 1;
                // Write FASTQ record
                writeln!(writer, "@{}", record.id)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&record.sequence))?;
                writeln!(writer, "+")?;
                writeln!(writer, "{}", String::from_utf8_lossy(&record.quality))?;
            }
        }

        let reads_failed = total_reads - reads_passed;
        let pass_rate = if total_reads > 0 {
            (reads_passed as f64 / total_reads as f64) * 100.0
        } else {
            0.0
        };

        Ok(QualityFilterStats {
            sample_name: sample_name.to_string(),
            total_reads,
            reads_passed,
            reads_failed,
            pass_rate,
            quality_threshold: self.min_quality,
        })
    }
}