//! Biometal QC Tools
//!
//! Fast QC tools for virome analysis using proven biometal primitives.
//!
//! This library provides shared functionality for:
//! - Quality statistics and base composition analysis
//! - Contamination screening
//! - VLP assessment metrics
//! - Multi-sample QC reporting

pub mod quality;
pub mod contamination;
pub mod vlp;
pub mod reporting;

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Common QC statistics structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcStats {
    pub sample_name: String,
    pub total_reads: u64,
    pub total_bases: u64,
    pub gc_content: f64,
    pub mean_quality: f64,
    pub q30_bases: f64,
    pub complexity_score: f64,
}

/// Read a FASTQ file and return basic metadata
pub fn get_file_info<P: AsRef<Path>>(path: P) -> Result<(String, u64)> {
    let filename = path
        .as_ref()
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string();

    let metadata = std::fs::metadata(&path)?;
    Ok((filename, metadata.len()))
}

/// Validate that biometal dependency is available
pub fn check_biometal_availability() -> bool {
    // This will be expanded when we integrate with biometal
    true
}