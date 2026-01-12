//! VLP (Virus-Like Particle) assessment module
//!
//! Composition-based VLP success metrics using biometal primitives:
//! - gc_content: For GC distribution analysis
//! - complexity: For sequence diversity assessment
//! - base_counting: For composition patterns

use anyhow::Result;
use biometal::io::{DataSource, FastqStream};
use biometal::operations::{complexity_score, gc_content};
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VlpReport {
    pub sample_name: String,
    pub total_reads: u64,
    pub gc_distribution_score: f64,
    pub complexity_diversity: f64,
    pub compositional_evenness: f64,
    pub vlp_success_score: f64,
}

/// VLP assessor using composition-based metrics
pub struct VlpAssessor {
    pub min_complexity: f64,
    pub optimal_gc_range: (f64, f64),
    pub min_length: usize,
}

impl Default for VlpAssessor {
    fn default() -> Self {
        Self {
            min_complexity: 0.7,
            optimal_gc_range: (0.35, 0.65), // Typical viral GC range
            min_length: 50,
        }
    }
}

impl VlpAssessor {
    pub fn new(min_complexity: f64, optimal_gc_range: (f64, f64), min_length: usize) -> Self {
        Self {
            min_complexity,
            optimal_gc_range,
            min_length,
        }
    }

    /// Assess VLP success using composition-based metrics
    /// Uses biometal gc_content and complexity primitives
    pub fn assess_vlp<P: AsRef<Path>>(&self, fastq_path: P) -> Result<VlpReport> {
        let sample_name = fastq_path
            .as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Initialize metrics collection
        let mut total_reads = 0u64;
        let mut gc_values = Vec::new();
        let mut complexity_values = Vec::new();
        let mut base_counts = [0u64; 4]; // A, T, G, C

        // Create biometal data source and stream
        let data_source = DataSource::from_path(&fastq_path);
        let fastq_stream = FastqStream::new(data_source)?;

        // Process reads using biometal streaming
        for record_result in fastq_stream {
            let record = record_result?;

            // Skip empty or short records
            if record.is_empty() || record.sequence.len() < self.min_length {
                continue;
            }

            total_reads += 1;

            // Use biometal gc_content primitive
            let gc_ratio = gc_content(&record.sequence);
            gc_values.push(gc_ratio);

            // Use biometal complexity primitive
            let complexity = complexity_score(&record.sequence);
            complexity_values.push(complexity);

            // Count individual bases for compositional evenness
            for &base in &record.sequence {
                match base {
                    b'A' | b'a' => base_counts[0] += 1,
                    b'T' | b't' => base_counts[1] += 1,
                    b'G' | b'g' => base_counts[2] += 1,
                    b'C' | b'c' => base_counts[3] += 1,
                    _ => {}, // Ignore ambiguous bases
                }
            }
        }

        // Calculate VLP success metrics
        let gc_distribution_score = self.calculate_gc_distribution_score(&gc_values);
        let complexity_diversity = if !complexity_values.is_empty() {
            complexity_values.iter().sum::<f64>() / complexity_values.len() as f64
        } else {
            0.0
        };
        let compositional_evenness = self.calculate_compositional_evenness(&base_counts);
        let vlp_success_score = self.calculate_success_score(
            gc_distribution_score,
            complexity_diversity,
            compositional_evenness,
        );

        let report = VlpReport {
            sample_name,
            total_reads,
            gc_distribution_score,
            complexity_diversity,
            compositional_evenness,
            vlp_success_score,
        };

        Ok(report)
    }

    /// Calculate GC distribution score based on diversity within optimal range
    fn calculate_gc_distribution_score(&self, gc_values: &[f64]) -> f64 {
        if gc_values.is_empty() {
            return 0.0;
        }

        // Calculate proportion of reads within optimal GC range
        let in_range_count = gc_values
            .iter()
            .filter(|&&gc| gc >= self.optimal_gc_range.0 && gc <= self.optimal_gc_range.1)
            .count();

        let in_range_proportion = in_range_count as f64 / gc_values.len() as f64;

        // Calculate GC diversity (standard deviation)
        let mean_gc = gc_values.iter().sum::<f64>() / gc_values.len() as f64;
        let variance = gc_values
            .iter()
            .map(|&gc| (gc - mean_gc).powi(2))
            .sum::<f64>()
            / gc_values.len() as f64;
        let std_dev = variance.sqrt();

        // Score combines range adherence with diversity (normalized std dev)
        (in_range_proportion * 0.7) + (std_dev.min(0.2) / 0.2 * 0.3)
    }

    /// Calculate compositional evenness using Shannon evenness index
    fn calculate_compositional_evenness(&self, base_counts: &[u64; 4]) -> f64 {
        let total_bases: u64 = base_counts.iter().sum();
        if total_bases == 0 {
            return 0.0;
        }

        // Calculate Shannon entropy
        let mut entropy = 0.0;
        for &count in base_counts {
            if count > 0 {
                let proportion = count as f64 / total_bases as f64;
                entropy -= proportion * proportion.ln();
            }
        }

        // Normalize by maximum possible entropy (ln(4) for 4 bases)
        entropy / 4.0_f64.ln()
    }

    /// Calculate overall VLP success score
    pub fn calculate_success_score(
        &self,
        gc_score: f64,
        complexity: f64,
        evenness: f64,
    ) -> f64 {
        // Weighted combination of metrics
        (gc_score * 0.3) + (complexity * 0.4) + (evenness * 0.3)
    }

    /// Determine if VLP preparation was successful
    pub fn is_vlp_successful(&self, report: &VlpReport) -> bool {
        report.vlp_success_score >= 0.7
            && report.complexity_diversity >= self.min_complexity
    }
}