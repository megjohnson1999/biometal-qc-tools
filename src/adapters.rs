//! Adapter trimming module using biometal primitives
//!
//! Uses proven biometal primitives:
//! - AdapterDetector: Built-in Illumina adapter detection with 8-15× NEON speedup
//! - find_patterns: Multi-pattern matching for adapter detection
//! - trim_start/trim_end: Fixed-position trimming based on adapter positions
//! - FastqStream: Streaming I/O for constant memory usage

use crate::QcStatsMarker;
use anyhow::Result;
use biometal::alignment::{MotifFinder, MotifPattern, MotifMatch};
use biometal::io::{DataSource, FastqStream, FastqWriter};
use biometal::operations::{trim_start, trim_end};
use biometal::FastqRecord;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AdapterStats {
    pub total_reads: usize,
    pub reads_with_adapters: usize,
    pub adapters_found: HashMap<String, usize>,
    pub total_bases_trimmed: usize,
    pub average_trim_length: f64,
}

impl Default for AdapterStats {
    fn default() -> Self {
        Self {
            total_reads: 0,
            reads_with_adapters: 0,
            adapters_found: HashMap::new(),
            total_bases_trimmed: 0,
            average_trim_length: 0.0,
        }
    }
}

impl QcStatsMarker for AdapterStats {}

/// Adapter trimmer using biometal primitives
pub struct AdapterTrimmer {
    pub min_adapter_length: usize,
    pub min_overlap: usize,
    pub trim_both_ends: bool,
}

impl Default for AdapterTrimmer {
    fn default() -> Self {
        Self {
            min_adapter_length: 8,   // Minimum adapter match length
            min_overlap: 5,          // Minimum overlap to consider for trimming
            trim_both_ends: true,    // Check both 5' and 3' ends
        }
    }
}

impl AdapterTrimmer {
    /// Create a new adapter trimmer with custom parameters
    pub fn new(min_adapter_length: usize, min_overlap: usize, trim_both_ends: bool) -> Self {
        Self {
            min_adapter_length,
            min_overlap,
            trim_both_ends,
        }
    }

    /// Process FASTQ file and trim adapters
    pub fn process_fastq<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: Option<P>,
    ) -> Result<AdapterStats> {
        // Create motif finder with Illumina adapters (same as AdapterDetector::new_illumina)
        let patterns = vec![
            MotifPattern::new("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "Illumina Universal"),
            MotifPattern::new("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "Illumina Small RNA 3'"),
            MotifPattern::new("TGGAATTCTCGGGTGCCAAGG", "Illumina Small RNA 5'"),
        ];
        let motif_finder = MotifFinder::new(patterns, 60); // High stringency like AdapterDetector

        let mut stats = AdapterStats::default();

        let input_path_ref = input_path.as_ref();
        let data_source = DataSource::from_path(input_path_ref);
        let fastq_stream = FastqStream::new(data_source)?;

        // Create streaming writer if output is requested (constant memory usage)
        let mut writer = if let Some(ref output_path) = output_path {
            Some(FastqWriter::create(output_path)?)
        } else {
            None
        };

        // Process records in streaming fashion
        for record_result in fastq_stream {
            let record = record_result?;
            stats.total_reads += 1;

            // Find adapters in this record
            let matches = motif_finder.find_in_sequence(
                &record.id,
                &record.sequence
            );

            let mut was_trimmed = false;
            let mut bases_trimmed_this_read = 0;
            let mut trimmed_record = record; // Use original record, clone only if needed

            if !matches.is_empty() {
                stats.reads_with_adapters += 1;

                // Process each adapter match
                for adapter_match in matches {
                    // Count adapter occurrences
                    *stats.adapters_found.entry(adapter_match.motif_name.clone())
                        .or_insert(0) += 1;

                    // Determine trim positions based on adapter location
                    let trim_pos = self.calculate_trim_position(&adapter_match, trimmed_record.sequence.len());

                    if let Some((trim_start_pos, trim_end_pos)) = trim_pos {
                        // Apply trimming based on position - now we'll modify the record
                        if trim_start_pos > 0 {
                            trimmed_record = trim_start(&trimmed_record, trim_start_pos)?;
                            bases_trimmed_this_read += trim_start_pos;
                            was_trimmed = true;
                        }

                        if trim_end_pos > 0 && trimmed_record.sequence.len() > trim_end_pos {
                            let new_length = trimmed_record.sequence.len() - trim_end_pos;
                            trimmed_record = trim_end(&trimmed_record, new_length)?;
                            bases_trimmed_this_read += trim_end_pos;
                            was_trimmed = true;
                        }
                    }
                }
            }

            if was_trimmed {
                stats.total_bases_trimmed += bases_trimmed_this_read;
            }

            // Write record immediately if output is requested (streaming)
            if let Some(ref mut w) = writer {
                w.write_record(&trimmed_record)?;
            }
        }

        // Calculate average trim length
        if stats.reads_with_adapters > 0 {
            stats.average_trim_length = stats.total_bases_trimmed as f64 / stats.reads_with_adapters as f64;
        }

        // Finalize output stream if opened
        if let Some(w) = writer {
            w.finish()?;
        }

        Ok(stats)
    }

    /// Calculate trim positions based on adapter match
    fn calculate_trim_position(&self, adapter_match: &MotifMatch, sequence_length: usize) -> Option<(usize, usize)> {
        let mut trim_start = 0;
        let mut trim_end = 0;

        // Adapter at the beginning (5' end)
        if adapter_match.position <= self.min_overlap {
            trim_start = adapter_match.position + adapter_match.length;
        }

        // Adapter at the end (3' end) - only check if trim_both_ends is true
        if self.trim_both_ends && adapter_match.position + adapter_match.length >= sequence_length - self.min_overlap {
            trim_end = sequence_length - adapter_match.position;
        }

        if trim_start > 0 || trim_end > 0 {
            Some((trim_start, trim_end))
        } else {
            None
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adapter_trimmer_creation() {
        let trimmer = AdapterTrimmer::default();
        assert_eq!(trimmer.min_adapter_length, 8);
        assert_eq!(trimmer.min_overlap, 5);
        assert!(trimmer.trim_both_ends);
    }

    #[test]
    fn test_adapter_trimmer_custom() {
        let trimmer = AdapterTrimmer::new(10, 3, false);
        assert_eq!(trimmer.min_adapter_length, 10);
        assert_eq!(trimmer.min_overlap, 3);
        assert!(!trimmer.trim_both_ends);
    }
}