//! Biometal PolyG Trimmer Tool
//!
//! Remove polyG tails from NovaSeq FASTQ reads caused by 2-channel chemistry artifacts.

use anyhow::Result;
use biometal::io::{DataSource, FastqStream};
use biometal_qc_tools::{get_file_info, PolyGStats};
use clap::{Arg, Command};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn main() -> Result<()> {
    let matches = Command::new("biometal-polyg-trim")
        .version("0.1.0")
        .about("Remove polyG tails from NovaSeq FASTQ reads")
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
                .help("Output trimmed FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("min_polyg_length")
                .long("min-polyg-length")
                .value_name("LENGTH")
                .help("Minimum consecutive Gs to trigger trimming")
                .default_value("10"),
        )
        .arg(
            Arg::new("min_read_length")
                .long("min-read-length")
                .value_name("LENGTH")
                .help("Minimum read length after trimming")
                .default_value("50"),
        )
        .arg(
            Arg::new("stats")
                .long("stats")
                .value_name("JSON")
                .help("Output statistics JSON file")
                .default_value("polyg_stats.json"),
        )
        .get_matches();

    // Parse arguments
    let input_file = PathBuf::from(matches.get_one::<String>("input").unwrap());
    let output_file = PathBuf::from(matches.get_one::<String>("output").unwrap());
    let min_polyg_length: usize = matches
        .get_one::<String>("min_polyg_length")
        .unwrap()
        .parse()?;
    let min_read_length: usize = matches
        .get_one::<String>("min_read_length")
        .unwrap()
        .parse()?;
    let stats_file = PathBuf::from(matches.get_one::<String>("stats").unwrap());

    println!("✂️  Biometal PolyG Trimmer");
    println!("Input: {}", input_file.display());
    println!("Output: {}", output_file.display());
    println!("Min polyG length: {}", min_polyg_length);
    println!("Min read length: {}", min_read_length);

    // Get sample name from input file
    let (sample_name, _) = get_file_info(&input_file)?;

    // Create trimmer and process
    let trimmer = PolyGTrimmer::new(min_polyg_length, min_read_length);
    let stats = trimmer.trim_reads(&input_file, &output_file, &sample_name)?;

    // Output statistics
    println!("📊 PolyG Trimming Results:");
    println!("  Total reads: {}", stats.total_reads);
    println!("  Reads trimmed: {} ({:.1}%)",
             stats.reads_trimmed,
             (stats.reads_trimmed as f64 / stats.total_reads as f64) * 100.0);
    println!("  Reads discarded: {} ({:.1}%)",
             stats.reads_discarded,
             (stats.reads_discarded as f64 / stats.total_reads as f64) * 100.0);
    println!("  Total bases removed: {}", stats.total_bases_removed);
    println!("  Average trim length: {:.1} bases", stats.average_trim_length);

    // Export statistics to JSON
    let json_content = serde_json::to_string_pretty(&stats)?;
    std::fs::write(&stats_file, json_content)?;
    println!("💾 Statistics saved to: {}", stats_file.display());

    Ok(())
}

/// PolyG trimmer implementation
pub struct PolyGTrimmer {
    pub min_polyg_length: usize,
    pub min_read_length: usize,
}

impl PolyGTrimmer {
    pub fn new(min_polyg_length: usize, min_read_length: usize) -> Self {
        Self {
            min_polyg_length,
            min_read_length,
        }
    }

    /// Trim polyG tails from FASTQ reads
    pub fn trim_reads(
        &self,
        input_path: &PathBuf,
        output_path: &PathBuf,
        sample_name: &str,
    ) -> Result<PolyGStats> {
        let mut total_reads = 0u64;
        let mut reads_trimmed = 0u64;
        let mut reads_discarded = 0u64;
        let mut total_bases_removed = 0u64;

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

            // Trim polyG tail
            let (trimmed_sequence, trimmed_quality, trim_length) =
                self.trim_polyg_tail(&record.sequence, &record.quality);

            if trim_length > 0 {
                reads_trimmed += 1;
                total_bases_removed += trim_length as u64;
            }

            // Check if read meets minimum length requirement
            if trimmed_sequence.len() >= self.min_read_length {
                // Write trimmed FASTQ record
                writeln!(writer, "@{}", record.id)?;
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_sequence))?;
                writeln!(writer, "+")?;
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_quality))?;
            } else {
                reads_discarded += 1;
            }
        }

        let average_trim_length = if reads_trimmed > 0 {
            total_bases_removed as f64 / reads_trimmed as f64
        } else {
            0.0
        };

        Ok(PolyGStats {
            sample_name: sample_name.to_string(),
            total_reads,
            reads_trimmed,
            reads_discarded,
            total_bases_removed,
            average_trim_length,
        })
    }

    /// Trim polyG tail from 3' end of read
    fn trim_polyg_tail(&self, sequence: &[u8], quality: &[u8]) -> (Vec<u8>, Vec<u8>, usize) {
        let seq_len = sequence.len();

        // Scan from 3' end for consecutive Gs
        let mut polyg_start = seq_len;
        let mut consecutive_gs = 0;

        for i in (0..seq_len).rev() {
            if sequence[i] == b'G' || sequence[i] == b'g' {
                consecutive_gs += 1;
                if consecutive_gs >= self.min_polyg_length {
                    polyg_start = i;
                }
            } else {
                // Reset if we encounter a non-G base
                consecutive_gs = 0;
            }
        }

        if polyg_start < seq_len {
            // Trim back to before the polyG tail
            let trim_length = seq_len - polyg_start;
            let trimmed_seq = sequence[0..polyg_start].to_vec();
            let trimmed_qual = quality[0..polyg_start].to_vec();
            (trimmed_seq, trimmed_qual, trim_length)
        } else {
            // No polyG tail found, return original
            (sequence.to_vec(), quality.to_vec(), 0)
        }
    }
}