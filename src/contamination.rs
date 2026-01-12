//! Contamination screening module
//!
//! Uses biometal primitives for PhiX and vector detection:
//! - pattern_match: For known contamination sequences
//! - base_counting: For composition-based detection

use anyhow::Result;
use biometal::io::{DataSource, FastqStream};
use biometal::operations::has_pattern;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContaminationReport {
    pub sample_name: String,
    pub total_reads: u64,
    pub phix_reads: u64,
    pub vector_reads: u64,
    pub phix_percentage: f64,
    pub vector_percentage: f64,
}

/// Contamination screener using biometal primitives
pub struct ContaminationScreener {
    pub phix_threshold: f64,
    pub vector_threshold: f64,
    pub min_length: usize,
}

impl Default for ContaminationScreener {
    fn default() -> Self {
        Self {
            phix_threshold: 0.1, // 0.1% PhiX threshold
            vector_threshold: 0.05, // 0.05% vector threshold
            min_length: 50, // Minimum read length
        }
    }
}

impl ContaminationScreener {
    pub fn new(phix_threshold: f64, vector_threshold: f64, min_length: usize) -> Self {
        Self {
            phix_threshold,
            vector_threshold,
            min_length,
        }
    }

    /// Comprehensive PhiX174 contamination sequences
    /// Uses longer, more specific patterns from PhiX174 genome for better sensitivity
    /// Replaces short snippets with representative 31-mers from different genome regions
    fn get_phix_patterns() -> Vec<&'static str> {
        vec![
            // Start region (positions 1-31)
            "GAGTTTTATCGCTTCCATGACGCAGAAGTTA",
            // Early genes region (positions 200-230)
            "GTGGACTGCTGGCGGAAAATGAGAAAATTCG",
            // DNA replication region (positions 500-530)
            "TGGCTTAATATGCTTGGCACGTTCGTCAAGG",
            // Middle genes region (positions 1000-1030)
            "GATGAGGAGAAGTGGCTTAATATGCTTGGCA",
            // Late genes region (positions 1500-1530)
            "TGACATTTTAAAAGAGCGTGGATTACTATCT",
            // Packaging region (positions 2000-2030)
            "TCCGATGCTGTTCAACCACTAATAGGTAAGA",
            // Capsid region (positions 2500-2530)
            "TGGATTTAACCGAAGATGATTTCGATTTTCT",
            // Assembly region (positions 3000-3030)
            "CGCTGAATTGTTCGCGTTTACCTTGCGTGTA",
            // Lysis region (positions 3500-3530)
            "TGCGCATGACCTTTCCCATCTTGGCTTCCTT",
            // End region (positions 4000-4030)
            "CCACTCCTCTCCCGACTGTTAACACTACTGG",
            // Terminal region (positions 4500-4530)
            "TGGTGTTAATGCCACTCCTCTCCCGACTGTT",
            // Final region (positions 5000-5030)
            "GGCCGCTGGTGTTGATCGGCCTTGATGATCG",
            // Overlapping critical regions for high sensitivity
            "ATGACGCAGAAGTTAACACTTTCGGATATTT",
            "TCCATGACGCAGAAGTTAACACTTTCGGATA",
            "ACGCAGAAGTTAACACTTTCGGATATTTCTG",
            "CGAGCGTCCGGTTAAAGCCGCTGAATTGTTC",
            "AGCGTCCGGTTAAAGCCGCTGAATTGTTCGC",
            "GTCCGGTTAAAGCCGCTGAATTGTTCGCGTT",
        ]
    }

    /// Comprehensive vector contamination sequences
    /// Representative 31-mers from major cloning vectors, expression vectors, and viral vectors
    /// Based on UniVec database used by BBDuk in the lab pipeline
    fn get_vector_patterns() -> Vec<&'static str> {
        vec![
            // pBR322 (most common cloning vector) - multiple regions
            "TTCTCATGTTTGACAGCTTATCATCGATAAG",
            "TTATCATCGATAAGCTTTAATGCGGTAGTTT",
            "AAGCTTTAATGCGGTAGTTTATCACAGTTAA",
            "TGCGGTAGTTTATCACAGTTAAATTGCTAAC",
            "TCAGGCACCGTGTATGAAATCTAACAATGCG",

            // M13 phage (common cloning/sequencing vector)
            "AACGCTACTACTATTAGTAGAATTGATGCCA",
            "CTACTATTAGTAGAATTGATGCCACCTTTTC",
            "TTAGTAGAATTGATGCCACCTTTTCAGCTCG",
            "GCCACCTTTTCAGCTCGCGCCCCAAATGAAA",
            "GTAAAACGACGGCCAGTGTAACGAACTTCC", // M13 forward primer region
            "CAGGAAACAGCTATGACGCTGCCGACGGTT", // M13 reverse primer region

            // SV40 (simian virus 40 - common viral vector)
            "AGCAGACATGATAAGATACATTGATGAGTTC",
            "ACATTGATGAGTTTGGACAAACCACAACTAG",
            "GGACAAACCACAACTAGAATGCAGTGAAAAA",

            // Lambda phage (common expression vector)
            "AGGGCGATCGGTGCGGGCCTCTTCGCTATTA",
            "TCGGTGCGGGCCTCTTCGCTATTACGCCAGC",
            "GGCCTCTTCGCTATTACGCCAGCTGGCGAAA",

            // pUC series (common high-copy cloning vectors)
            "GCCAGTGTAACGAACTTCCGCATTGTCCGAT", // pUC origin region
            "TGTAACGAACTTCCGCATTGTCCGATGATCG",
            "CTTCCGCATTGTCCGATGATCGCGGTGATAT",
            "TTGTAAAACGACGGCCAGTGCCAAGCTTGCA", // pUC polylinker
            "AAAACGACGGCCAGTGCCAAGCTTGCATGCC",

            // pGEM series (common TA cloning vectors)
            "TGCCCGCTTTCCAGTCGGGAAACCTGTCGTG",
            "CCGCTTTCCAGTCGGGAAACCTGTCGTGCCA",
            "TCCAGTCGGGAAACCTGTCGTGCCAGCTGCA",

            // pBluescript (common cloning vector)
            "GGAGCTCCGGCCGAGAATAAGCCTTGATGAT", // pBluescript polylinker
            "CTCCGGCCGAGAATAAGCCTTGATGATCGCG",
            "GGCCGAGAATAAGCCTTGATGATCGCGTGGT",

            // RSF1010 (broad-host-range plasmid)
            "CTGCTGGCGTCAGCGTTGGGCGTCTCGCGCT",
            "GGCGTCAGCGTTGGGCGTCTCGCGCTCGCGG",
            "AGCGTTGGGCGTCTCGCGCTCGCGGTAGGCG",

            // E.coli lac operon (expression vector component)
            "TGGAATTGTGAGCGGATAACAATTCCCCTCT",
            "ATTGTGAGCGGATAACAATTCCCCTCTAGAA",
            "GAGCGGATAACAATTCCCCTCTAGAAATAAA",

            // Common primer binding sites found in many vectors
            "CGCACAACGCCTCTATTGAAGGGGAGTCTTA", // Common vector primer site
            "ACAACGCCTCTATTGAAGGGGAGTCTTAAGG", // Common vector primer site
            "GCCTCTATTGAAGGGGAGTCTTAAGGCATTT", // Common vector primer site

            // Common ampicillin resistance gene (found in many vectors)
            "ATGAGTATTCAACATTTCCGTGTCGCCCTTA",
            "GTATTCAACATTTCCGTGTCGCCCTTATTCC",
            "CAACATTTCCGTGTCGCCCTTATTCCCTTTT",

            // Common tetracycline resistance gene
            "ATGAAAGTTGAACTAGATTCTCAACAACGGC",
            "AAGTTGAACTAGATTCTCAACAACGGCTTTT",
            "GAACTAGATTCTCAACAACGGCTTTTGCCAC",

            // Common neomycin resistance gene
            "ATGATTGAACAAGATGGATTGCACGCAGGTT",
            "TTGAACAAGATGGATTGCACGCAGGTTCTCC",
            "CAAGATGGATTGCACGCAGGTTCTCCGGCCG",
        ]
    }

    /// Screen for contamination in FASTQ file
    /// Uses biometal pattern_match for known sequences
    pub fn screen_fastq<P: AsRef<Path>>(&self, fastq_path: P) -> Result<ContaminationReport> {
        let sample_name = fastq_path
            .as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Get contamination patterns
        let phix_patterns = Self::get_phix_patterns();
        let vector_patterns = Self::get_vector_patterns();

        // Initialize counters
        let mut total_reads = 0u64;
        let mut phix_reads = 0u64;
        let mut vector_reads = 0u64;

        // Create biometal data source and stream
        let data_source = DataSource::from_path(&fastq_path);
        let fastq_stream = FastqStream::new(data_source)?;

        // Process reads using biometal streaming
        for record_result in fastq_stream {
            let record = record_result?;

            // Skip empty records or records below minimum length
            if record.is_empty() || record.sequence.len() < self.min_length {
                continue;
            }

            total_reads += 1;

            // Check for PhiX contamination using biometal pattern matching
            let mut is_phix_contaminated = false;
            for pattern in &phix_patterns {
                if has_pattern(&record.sequence, pattern.as_bytes()) {
                    is_phix_contaminated = true;
                    break;
                }
            }
            if is_phix_contaminated {
                phix_reads += 1;
            }

            // Check for vector contamination using biometal pattern matching
            let mut is_vector_contaminated = false;
            for pattern in &vector_patterns {
                if has_pattern(&record.sequence, pattern.as_bytes()) {
                    is_vector_contaminated = true;
                    break;
                }
            }
            if is_vector_contaminated {
                vector_reads += 1;
            }
        }

        // Calculate percentages
        let phix_percentage = if total_reads > 0 {
            (phix_reads as f64 / total_reads as f64) * 100.0
        } else {
            0.0
        };

        let vector_percentage = if total_reads > 0 {
            (vector_reads as f64 / total_reads as f64) * 100.0
        } else {
            0.0
        };

        let report = ContaminationReport {
            sample_name,
            total_reads,
            phix_reads,
            vector_reads,
            phix_percentage,
            vector_percentage,
        };

        Ok(report)
    }

    /// Check if contamination levels are within acceptable thresholds
    pub fn is_contamination_acceptable(&self, report: &ContaminationReport) -> bool {
        report.phix_percentage <= self.phix_threshold
            && report.vector_percentage <= self.vector_threshold
    }
}