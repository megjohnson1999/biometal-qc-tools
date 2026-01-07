//! Quality statistics module using biometal primitives
//!
//! Uses proven biometal primitives:
//! - base_counting: For base composition analysis
//! - gc_content: For GC content calculation
//! - quality_filter: For quality score analysis
//! - complexity: For sequence complexity assessment

use crate::QcStats;
use anyhow::Result;
use biometal::io::{DataSource, FastqStream};
use biometal::operations::{complexity_score, gc_content, mean_quality};
use std::path::Path;

/// Quality statistics calculator using biometal primitives
pub struct QualityAnalyzer {
    pub min_quality: u8,
    pub min_length: usize,
}

impl Default for QualityAnalyzer {
    fn default() -> Self {
        Self {
            min_quality: 20,
            min_length: 50,
        }
    }
}

impl QualityAnalyzer {
    pub fn new(min_quality: u8, min_length: usize) -> Self {
        Self {
            min_quality,
            min_length,
        }
    }

    /// Analyze quality statistics from a FASTQ file
    /// Uses biometal primitives: base_counting, gc_content, quality_filter, complexity
    pub fn analyze_fastq<P: AsRef<Path>>(&self, fastq_path: P) -> Result<QcStats> {
        let sample_name = fastq_path
            .as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Initialize counters
        let mut total_reads = 0u64;
        let mut total_bases = 0u64;
        let mut gc_count = 0u64;
        let mut quality_sum = 0f64;
        let mut q30_bases_count = 0u64;
        let mut complexity_sum = 0f64;

        // Create biometal data source and stream
        let data_source = DataSource::from_path(&fastq_path);
        let fastq_stream = FastqStream::new(data_source)?;

        // Process records using biometal streaming
        for record_result in fastq_stream {
            let record = record_result?;

            // Skip empty records
            if record.is_empty() || record.sequence.len() < self.min_length {
                continue;
            }

            total_reads += 1;
            total_bases += record.sequence.len() as u64;

            // Use biometal gc_content primitive
            let gc_content_ratio = gc_content(&record.sequence);
            gc_count += (gc_content_ratio * record.sequence.len() as f64) as u64;

            // Use biometal mean_quality primitive
            let record_mean_quality = mean_quality(&record.quality);
            quality_sum += record_mean_quality;

            // Count Q30 bases (quality >= 30, which is 63 in Phred+33)
            let q30_count = record.quality.iter().filter(|&&q| q >= 63).count();
            q30_bases_count += q30_count as u64;

            // Use biometal complexity primitive
            let record_complexity = complexity_score(&record.sequence);
            complexity_sum += record_complexity;
        }

        // Calculate final statistics
        let gc_content_percent = if total_bases > 0 {
            (gc_count as f64 / total_bases as f64) * 100.0
        } else {
            0.0
        };

        let mean_quality_score = if total_reads > 0 {
            quality_sum / total_reads as f64
        } else {
            0.0
        };

        let q30_percentage = if total_bases > 0 {
            (q30_bases_count as f64 / total_bases as f64) * 100.0
        } else {
            0.0
        };

        let avg_complexity = if total_reads > 0 {
            complexity_sum / total_reads as f64
        } else {
            0.0
        };

        let stats = QcStats {
            sample_name,
            total_reads,
            total_bases,
            gc_content: gc_content_percent,
            mean_quality: mean_quality_score,
            q30_bases: q30_percentage,
            complexity_score: avg_complexity,
        };

        Ok(stats)
    }

    /// Calculate quality distribution metrics
    pub fn quality_distribution(&self, qualities: &[u8]) -> QualityDistribution {
        let total = qualities.len();
        // Convert Phred+33 to actual quality scores for thresholds
        let q30_count = qualities.iter().filter(|&&q| q >= 63).count(); // Q30 = 30 + 33
        let q20_count = qualities.iter().filter(|&&q| q >= 53).count(); // Q20 = 20 + 33

        QualityDistribution {
            q30_percent: (q30_count as f64 / total as f64) * 100.0,
            q20_percent: (q20_count as f64 / total as f64) * 100.0,
            mean_quality: mean_quality(qualities),
        }
    }
}

#[derive(Debug, Clone)]
pub struct QualityDistribution {
    pub q30_percent: f64,
    pub q20_percent: f64,
    pub mean_quality: f64,
}