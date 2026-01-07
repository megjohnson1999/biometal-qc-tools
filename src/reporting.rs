//! Multi-sample QC reporting module
//!
//! Aggregates and reports QC metrics across multiple samples

use crate::{contamination::ContaminationReport, QcStats, vlp::VlpReport};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleQcReport {
    pub quality_stats: QcStats,
    pub contamination_report: ContaminationReport,
    pub vlp_report: VlpReport,
    pub overall_pass: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiSampleReport {
    pub samples: Vec<SampleQcReport>,
    pub summary: QcSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcSummary {
    pub total_samples: usize,
    pub passed_samples: usize,
    pub failed_samples: usize,
    pub pass_rate: f64,
    pub average_quality: f64,
    pub average_gc_content: f64,
}

/// QC reporter for multi-sample analysis
pub struct QcReporter {
    pub quality_threshold: f64,
    pub contamination_threshold: f64,
}

impl Default for QcReporter {
    fn default() -> Self {
        Self {
            quality_threshold: 25.0,
            contamination_threshold: 0.1,
        }
    }
}

impl QcReporter {
    pub fn new(quality_threshold: f64, contamination_threshold: f64) -> Self {
        Self {
            quality_threshold,
            contamination_threshold,
        }
    }

    /// Generate comprehensive QC report for multiple samples
    pub fn generate_report(&self, sample_reports: Vec<SampleQcReport>) -> MultiSampleReport {
        let total_samples = sample_reports.len();
        let passed_samples = sample_reports.iter().filter(|r| r.overall_pass).count();
        let failed_samples = total_samples - passed_samples;

        let average_quality = if !sample_reports.is_empty() {
            sample_reports
                .iter()
                .map(|r| r.quality_stats.mean_quality)
                .sum::<f64>()
                / total_samples as f64
        } else {
            0.0
        };

        let average_gc_content = if !sample_reports.is_empty() {
            sample_reports
                .iter()
                .map(|r| r.quality_stats.gc_content)
                .sum::<f64>()
                / total_samples as f64
        } else {
            0.0
        };

        let summary = QcSummary {
            total_samples,
            passed_samples,
            failed_samples,
            pass_rate: (passed_samples as f64 / total_samples as f64) * 100.0,
            average_quality,
            average_gc_content,
        };

        MultiSampleReport {
            samples: sample_reports,
            summary,
        }
    }

    /// Export report to JSON
    pub fn export_json<P: AsRef<Path>>(&self, report: &MultiSampleReport, path: P) -> Result<()> {
        let json_content = serde_json::to_string_pretty(report)?;
        std::fs::write(path, json_content)?;
        Ok(())
    }

    /// Determine if sample passes overall QC
    pub fn evaluate_sample(&self, sample: &SampleQcReport) -> bool {
        sample.quality_stats.mean_quality >= self.quality_threshold
            && sample.contamination_report.phix_percentage <= self.contamination_threshold
            && sample.vlp_report.vlp_success_score >= 0.7
    }
}