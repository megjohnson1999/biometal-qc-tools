//! Biometal Optical Duplicate Detection Tool
//!
//! Fast optical duplicate detection using biometal spatial primitives
//! Uses NEON-optimized coordinate processing to replace clumpify optical deduplication

use anyhow::Result;
use biometal::operations::spatial::{
    parse_illumina_coordinates, find_optical_duplicates,
    IlluminaCoordinate
};
use biometal::{FastqStream, FastqWriter, FastqRecord};
use clap::{Arg, Command};
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Serialize, Deserialize)]
struct OpticalDedupStats {
    total_reads: u64,
    optical_duplicates_found: u64,
    unique_reads_kept: u64,
    duplicate_groups: u64,
    average_group_size: f64,
    distance_threshold: f64,
    processing_time_seconds: f64,
}

fn main() -> Result<()> {
    let matches = Command::new("biometal-optical-dedup")
        .version("0.1.0")
        .about("Fast optical duplicate detection using biometal spatial primitives")
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
                .help("Output deduplicated FASTQ file")
                .required(true),
        )
        .arg(
            Arg::new("stats")
                .short('s')
                .long("stats")
                .value_name("JSON")
                .help("Output optical deduplication statistics (JSON)")
                .default_value("optical_dedup_stats.json"),
        )
        .arg(
            Arg::new("threshold")
                .short('t')
                .long("threshold")
                .value_name("PIXELS")
                .help("Distance threshold for optical duplicates (pixels)")
                .default_value("10.0"),
        )
        .arg(
            Arg::new("keep_best")
                .long("keep-best-quality")
                .help("Keep highest quality read from each duplicate group")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    // Parse arguments
    let input_path: PathBuf = matches.get_one::<String>("input").unwrap().into();
    let output_path: PathBuf = matches.get_one::<String>("output").unwrap().into();
    let stats_path: PathBuf = matches.get_one::<String>("stats").unwrap().into();
    let threshold: f64 = matches.get_one::<String>("threshold").unwrap().parse()?;
    let keep_best_quality = matches.get_flag("keep_best");

    println!("🔬 Biometal Optical Duplicate Detection");
    println!("=====================================");
    println!("Input: {}", input_path.display());
    println!("Output: {}", output_path.display());
    println!("Distance threshold: {} pixels", threshold);
    println!("Quality selection: {}", if keep_best_quality { "Best quality" } else { "First occurrence" });

    let start_time = std::time::Instant::now();

    // Process optical duplicates
    let stats = process_optical_duplicates(
        &input_path,
        &output_path,
        threshold,
        keep_best_quality,
    )?;

    let processing_time = start_time.elapsed().as_secs_f64();
    let final_stats = OpticalDedupStats {
        processing_time_seconds: processing_time,
        ..stats
    };

    // Write statistics
    let stats_json = serde_json::to_string_pretty(&final_stats)?;
    std::fs::write(&stats_path, stats_json)?;

    println!("\n✅ Optical Deduplication Complete");
    println!("Total reads processed: {}", final_stats.total_reads);
    println!("Optical duplicates found: {}", final_stats.optical_duplicates_found);
    println!("Unique reads kept: {}", final_stats.unique_reads_kept);
    println!("Duplicate groups: {}", final_stats.duplicate_groups);
    println!("Average group size: {:.2}", final_stats.average_group_size);
    println!("Processing time: {:.2}s", processing_time);
    println!("Statistics written to: {}", stats_path.display());

    Ok(())
}

fn process_optical_duplicates(
    input_path: &PathBuf,
    output_path: &PathBuf,
    threshold: f64,
    keep_best_quality: bool,
) -> Result<OpticalDedupStats> {
    // Step 1: Read all records and extract coordinates
    let mut fastq_stream = FastqStream::from_path(input_path)?;
    let mut records_with_coords: Vec<(FastqRecord, IlluminaCoordinate)> = Vec::new();
    let mut total_reads = 0;
    let mut parse_errors = 0;

    println!("📊 Reading FASTQ and parsing coordinates...");

    for record_result in fastq_stream {
        let record = record_result?;
        total_reads += 1;

        match parse_illumina_coordinates(&record.id) {
            Ok(coord) => {
                records_with_coords.push((record, coord));
            }
            Err(_) => {
                parse_errors += 1;
                // Keep records that can't be parsed (non-Illumina format)
                records_with_coords.push((record, create_default_coordinate()));
            }
        }
    }

    println!("   - Total reads: {}", total_reads);
    println!("   - Parseable coordinates: {}", records_with_coords.len() - parse_errors);
    println!("   - Parse errors: {}", parse_errors);

    // Step 2: Extract coordinates for optical duplicate detection
    let coordinates: Vec<IlluminaCoordinate> = records_with_coords
        .iter()
        .map(|(_, coord)| coord.clone())
        .collect();

    println!("🔍 Finding optical duplicates...");

    // Step 3: Find optical duplicate groups using biometal spatial primitives
    let duplicate_groups = find_optical_duplicates(
        coordinates.into_iter(),
        threshold
    )?;

    println!("   - Found {} duplicate groups", duplicate_groups.len());

    // Step 4: Create duplicate index mapping
    let mut duplicate_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut representatives: HashMap<usize, usize> = HashMap::new();

    for group in &duplicate_groups {
        if group.len() > 1 {
            // Select representative (first or best quality)
            let representative = if keep_best_quality {
                select_best_quality_read(group, &records_with_coords)?
            } else {
                group[0] // First occurrence
            };

            representatives.insert(representative, group.len());

            // Mark all others as duplicates
            for &idx in group {
                if idx != representative {
                    duplicate_indices.insert(idx);
                }
            }
        }
    }

    println!("   - Optical duplicates to remove: {}", duplicate_indices.len());

    // Step 5: Write filtered output
    println!("📝 Writing deduplicated output...");

    let mut writer = FastqWriter::create(output_path)?;
    let mut unique_reads_kept = 0;

    for (i, (record, _)) in records_with_coords.iter().enumerate() {
        if !duplicate_indices.contains(&i) {
            writer.write_record(record)?;
            unique_reads_kept += 1;
        }
    }

    let average_group_size = if duplicate_groups.is_empty() {
        0.0
    } else {
        duplicate_groups.iter().map(|g| g.len()).sum::<usize>() as f64 / duplicate_groups.len() as f64
    };

    Ok(OpticalDedupStats {
        total_reads,
        optical_duplicates_found: duplicate_indices.len() as u64,
        unique_reads_kept,
        duplicate_groups: duplicate_groups.len() as u64,
        average_group_size,
        distance_threshold: threshold,
        processing_time_seconds: 0.0, // Will be set by caller
    })
}

fn select_best_quality_read(
    group: &[usize],
    records_with_coords: &[(FastqRecord, IlluminaCoordinate)]
) -> Result<usize> {
    let mut best_idx = group[0];
    let mut best_quality = calculate_mean_quality(&records_with_coords[group[0]].0)?;

    for &idx in &group[1..] {
        let quality = calculate_mean_quality(&records_with_coords[idx].0)?;
        if quality > best_quality {
            best_quality = quality;
            best_idx = idx;
        }
    }

    Ok(best_idx)
}

fn calculate_mean_quality(record: &FastqRecord) -> Result<f64> {
    let qualities = &record.quality;
    if qualities.is_empty() {
        return Ok(0.0);
    }

    let sum: u32 = qualities.iter().map(|&q| (q - 33) as u32).sum();
    Ok(sum as f64 / qualities.len() as f64)
}

// Helper function to create default coordinate for unparseable reads
fn create_default_coordinate() -> IlluminaCoordinate {
    IlluminaCoordinate {
        instrument: "UNKNOWN".to_string(),
        run_id: 0,
        flowcell: "UNKNOWN".to_string(),
        lane: 0,
        tile: 0,
        x: 0,
        y: 0,
        read: 1,
        filtered: false,
        control: 0,
        index: String::new(),
    }
}