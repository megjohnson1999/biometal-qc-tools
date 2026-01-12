//! Primer removal module using biometal primitives
//!
//! Replicates BBDuk's primer B removal functionality:
//! - Uses 24 Primer B variants from lab-virome-QC pipeline
//! - Implements k-mer matching (k=16, mink=9) using multiple pattern variants
//! - Two-step process: forward primers (5' trimming) + reverse complement primers (3' trimming)
//! - Uses biometal MotifFinder for pattern detection and trimming operations

use crate::QcStatsMarker;
use anyhow::Result;
use biometal::alignment::{MotifFinder, MotifPattern, MotifMatch};
use biometal::io::{DataSource, FastqStream};
use biometal::operations::{trim_start, trim_end};
use biometal::FastqRecord;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PrimerRemovalStats {
    pub sample_name: String,
    pub total_reads: usize,
    pub reads_with_forward_primers: usize,
    pub reads_with_rc_primers: usize,
    pub forward_primers_found: HashMap<String, usize>,
    pub rc_primers_found: HashMap<String, usize>,
    pub total_bases_trimmed: usize,
    pub contamination_level: f64, // Percentage of reads with unexpected primer variants
}

impl Default for PrimerRemovalStats {
    fn default() -> Self {
        Self {
            sample_name: String::new(),
            total_reads: 0,
            reads_with_forward_primers: 0,
            reads_with_rc_primers: 0,
            forward_primers_found: HashMap::new(),
            rc_primers_found: HashMap::new(),
            total_bases_trimmed: 0,
            contamination_level: 0.0,
        }
    }
}

impl QcStatsMarker for PrimerRemovalStats {}

/// Primer remover using biometal primitives
/// Replicates BBDuk's two-step primer B removal process
pub struct PrimerRemover {
    pub min_match_length: usize,  // Minimum k-mer size (BBDuk's mink=9)
    pub max_match_length: usize,  // Maximum k-mer size (BBDuk's k=16)
    pub contamination_threshold: f64, // Threshold for flagging cross-contamination
}

impl Default for PrimerRemover {
    fn default() -> Self {
        Self {
            min_match_length: 9,   // BBDuk mink=9
            max_match_length: 16,  // BBDuk k=16
            contamination_threshold: 0.05, // 5% contamination threshold
        }
    }
}

impl PrimerRemover {
    /// Create a new primer remover with custom parameters
    pub fn new(min_match_length: usize, max_match_length: usize, contamination_threshold: f64) -> Self {
        Self {
            min_match_length,
            max_match_length,
            contamination_threshold,
        }
    }

    /// Get 24 Primer B forward sequences from lab-virome-QC pipeline
    fn get_primer_b_forward_sequences() -> Vec<(&'static str, &'static str)> {
        vec![
            ("3GB-1", "TACCGTAGAGCTGCTA"),
            ("3GB-2", "ATAGAGCCTACTGTCG"),
            ("3GB-3", "GGGCCTTTAAGATCAC"),
            ("3GB-4", "CGTAGGAACGTCTCTA"),
            ("3GB-5", "TGAGCAGTGCATCATC"),
            ("3GB-6", "GCAACCTGGTCATGAT"),
            ("3GB-7", "ATGGAAGGTCATCTCC"),
            ("3GB-8", "CATCGGGCCATAATGT"),
            ("3GB-9", "CCTATTCATAGCGGGA"),
            ("3GB-10", "ACTTGATCGTCAACGG"),
            ("3GB-11", "CTCATCTGAGACGGAT"),
            ("3GB-12", "ATCTGCGCATGGATCA"),
            ("3GB-13", "GTAGTCAACTCCTGGA"),
            ("3GB-14", "AGGACGCTATGACTCT"),
            ("3GB-15", "CCATAGTGGCTGATCA"),
            ("3GB-16", "GGTTCACTTGAGACAC"),
            ("3GB-17", "ACGGTGCATCATTGAC"),
            ("3GB-18", "CATAGTGAGTCTGCCA"),
            ("3GB-19", "CTACATGCAGGATCTG"),
            ("3GB-20", "TACAACTGGATCGGTC"),
            ("3GB-21", "TCCGTAATCTGCAGAG"),
            ("3GB-22", "GTCTGTTCCAAACAGG"),
            ("3GB-23", "AGTCGCAGAGCTTCTA"),
            ("3GB-24", "TCACCAGGACATGTGT"),
        ]
    }

    /// Get 24 Primer B reverse complement sequences from lab-virome-QC pipeline
    fn get_primer_b_rc_sequences() -> Vec<(&'static str, &'static str)> {
        vec![
            ("3GB-1", "TAGCAGCTCTACGGTA"),
            ("3GB-2", "CGACAGTAGGCTCTAT"),
            ("3GB-3", "GTGATCTTAAAGGCCC"),
            ("3GB-4", "TAGAGACGTTCCTACG"),
            ("3GB-5", "GATGATGCACTGCTCA"),
            ("3GB-6", "ATCATGACCAGGTTGC"),
            ("3GB-7", "GGAGATGACCTTCCAT"),
            ("3GB-8", "ACATTATGGCCCGATG"),
            ("3GB-9", "TCCCGCTATGAATAGG"),
            ("3GB-10", "CCGTTGACGATCAAGT"),
            ("3GB-11", "ATCCGTCTCAGATGAG"),
            ("3GB-12", "TGATCCATGCGCAGAT"),
            ("3GB-13", "TCCAGGAGTTGACTAC"),
            ("3GB-14", "AGAGTCATAGCGTCCT"),
            ("3GB-15", "TGATCAGCCACTATGG"),
            ("3GB-16", "GTGTCTCAAGTGAACC"),
            ("3GB-17", "GTCAATGATGCACCGT"),
            ("3GB-18", "TGGCAGACTCACTATG"),
            ("3GB-19", "CAGATCCTGCATGTAG"),
            ("3GB-20", "GACCGATCCAGTTGTA"),
            ("3GB-21", "CTCTGCAGATTACGGA"),
            ("3GB-22", "CCTGTTTGGAACAGAC"),
            ("3GB-23", "TAGAAGCTCTGCGACT"),
            ("3GB-24", "ACACATGTCCTGGTGA"),
        ]
    }

    /// Generate k-mer variants for BBDuk-style matching (k=16 down to mink=9)
    /// This replicates BBDuk's ability to match partial primers
    fn generate_kmer_variants(sequence: &str, primer_id: &str, min_k: usize, max_k: usize) -> Vec<MotifPattern> {
        let mut patterns = Vec::new();

        // Generate k-mers from max_k down to min_k (16->15->14->...->9)
        for k in (min_k..=max_k).rev() {
            if sequence.len() >= k {
                // For forward primers: take k-mer from start (5' end)
                // For RC primers: take k-mer from end (3' end)
                let kmer = if k == sequence.len() {
                    sequence.to_string()
                } else {
                    sequence[0..k].to_string() // Always take from start for both forward and RC
                };

                let pattern_name = format!("{}_k{}", primer_id, k);
                patterns.push(MotifPattern::new(&kmer, &pattern_name));
            }
        }

        patterns
    }

    /// Create MotifFinder for forward primer detection (5' end matching)
    fn create_forward_primer_finder(&self) -> MotifFinder {
        let mut patterns = Vec::new();

        for (primer_id, sequence) in Self::get_primer_b_forward_sequences() {
            let kmer_patterns = Self::generate_kmer_variants(
                sequence,
                primer_id,
                self.min_match_length,
                self.max_match_length
            );
            patterns.extend(kmer_patterns);
        }

        // Use moderate threshold for primer detection (lower than adapters)
        MotifFinder::new(patterns, 30) // Moderate threshold for primer detection
    }

    /// Create MotifFinder for reverse complement primer detection (3' end matching)
    fn create_rc_primer_finder(&self) -> MotifFinder {
        let mut patterns = Vec::new();

        for (primer_id, sequence) in Self::get_primer_b_rc_sequences() {
            let kmer_patterns = Self::generate_kmer_variants(
                sequence,
                primer_id,
                self.min_match_length,
                self.max_match_length
            );
            patterns.extend(kmer_patterns);
        }

        // Use moderate threshold for primer detection (lower than adapters)
        MotifFinder::new(patterns, 30) // Moderate threshold for primer detection
    }

    /// Process FASTQ file and remove primers (two-step process like BBDuk)
    pub fn process_fastq<P: AsRef<Path>>(
        &self,
        input_path: P,
        output_path: Option<P>,
    ) -> Result<PrimerRemovalStats> {
        let sample_name = input_path
            .as_ref()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let mut stats = PrimerRemovalStats::default();
        stats.sample_name = sample_name;

        // Create primer finders
        let forward_finder = self.create_forward_primer_finder();
        let rc_finder = self.create_rc_primer_finder();

        let input_path_ref = input_path.as_ref();
        let data_source = DataSource::from_path(input_path_ref);
        let fastq_stream = FastqStream::new(data_source)?;

        let mut processed_records = Vec::new();

        // Process records in streaming fashion
        for record_result in fastq_stream {
            let record = record_result?;
            stats.total_reads += 1;

            // Step 1: Remove forward primers (5' end trimming like BBDuk ktrim="l")
            let mut current_record = record.clone();
            let forward_matches = forward_finder.find_in_sequence(&record.id, &record.sequence);

            if !forward_matches.is_empty() {
                stats.reads_with_forward_primers += 1;

                // Find the longest match (prioritize longer k-mers like BBDuk)
                if let Some(best_match) = self.find_best_forward_match(&forward_matches) {
                    // Extract primer ID from pattern name (e.g., "3GB-1_k16" -> "3GB-1")
                    let primer_id = best_match.motif_name.split('_').next().unwrap_or("unknown").to_string();
                    *stats.forward_primers_found.entry(primer_id).or_insert(0) += 1;

                    // Trim from 5' end (start of sequence)
                    let trim_length = best_match.position + best_match.length;
                    if trim_length < current_record.sequence.len() {
                        current_record = trim_start(&current_record, trim_length)?;
                        stats.total_bases_trimmed += trim_length;
                    }
                }
            }

            // Step 2: Remove reverse complement primers (3' end trimming like BBDuk ktrim="r")
            let rc_matches = rc_finder.find_in_sequence(&current_record.id, &current_record.sequence);

            if !rc_matches.is_empty() {
                stats.reads_with_rc_primers += 1;

                // Find the best match at 3' end
                if let Some(best_match) = self.find_best_rc_match(&rc_matches, current_record.sequence.len()) {
                    // Extract primer ID from pattern name
                    let primer_id = best_match.motif_name.split('_').next().unwrap_or("unknown").to_string();
                    *stats.rc_primers_found.entry(primer_id).or_insert(0) += 1;

                    // Trim from 3' end (end of sequence)
                    let new_length = best_match.position;
                    if new_length < current_record.sequence.len() && new_length > 0 {
                        let original_length = current_record.sequence.len();
                        current_record = trim_end(&current_record, new_length)?;
                        stats.total_bases_trimmed += original_length.saturating_sub(new_length);
                    }
                }
            }

            processed_records.push(current_record);
        }

        // Calculate contamination level (cross-contamination detection)
        stats.contamination_level = self.calculate_contamination_level(&stats);

        // Write output if requested
        if let Some(output_path) = output_path {
            self.write_trimmed_fastq(&processed_records, output_path)?;
        }

        Ok(stats)
    }

    /// Find the best forward primer match (longest k-mer at 5' end)
    fn find_best_forward_match<'a>(&self, matches: &'a [MotifMatch]) -> Option<&'a MotifMatch> {
        // Prioritize matches at the very beginning (position 0 or near it)
        // Among those, choose the longest k-mer
        matches.iter()
            .filter(|m| m.position <= 2) // Allow slight offset for sequencing errors
            .max_by_key(|m| {
                // Extract k-mer size from pattern name (e.g., "3GB-1_k16" -> 16)
                m.motif_name.split('_')
                    .nth(1)
                    .and_then(|s| s.strip_prefix('k'))
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(0)
            })
    }

    /// Find the best reverse complement primer match (longest k-mer at 3' end)
    fn find_best_rc_match<'a>(&self, matches: &'a [MotifMatch], sequence_length: usize) -> Option<&'a MotifMatch> {
        // Prioritize matches near the end of the sequence
        // Among those, choose the longest k-mer
        matches.iter()
            .filter(|m| {
                let distance_from_end = sequence_length.saturating_sub(m.position + m.length);
                distance_from_end <= 2 // Allow slight offset for sequencing errors
            })
            .max_by_key(|m| {
                // Extract k-mer size from pattern name
                m.motif_name.split('_')
                    .nth(1)
                    .and_then(|s| s.strip_prefix('k'))
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(0)
            })
    }

    /// Calculate contamination level based on primer diversity
    /// High diversity indicates cross-contamination between samples
    fn calculate_contamination_level(&self, stats: &PrimerRemovalStats) -> f64 {
        let total_primer_reads = stats.reads_with_forward_primers + stats.reads_with_rc_primers;
        if total_primer_reads == 0 {
            return 0.0;
        }

        // Calculate diversity (number of different primer variants found)
        let mut all_primers = HashMap::new();
        for (primer, count) in &stats.forward_primers_found {
            *all_primers.entry(primer.clone()).or_insert(0) += count;
        }
        for (primer, count) in &stats.rc_primers_found {
            *all_primers.entry(primer.clone()).or_insert(0) += count;
        }

        // If more than 1 primer variant is found, it might indicate contamination
        if all_primers.len() > 1 {
            // Calculate percentage of reads with non-dominant primer variants
            let max_count = all_primers.values().max().unwrap_or(&0);
            let minority_count: usize = all_primers.values().filter(|&&c| c != *max_count).sum();
            (minority_count as f64 / total_primer_reads as f64) * 100.0
        } else {
            0.0
        }
    }

    /// Write trimmed FASTQ records to file
    fn write_trimmed_fastq<P: AsRef<Path>>(
        &self,
        records: &[FastqRecord],
        output_path: P,
    ) -> Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};

        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        for record in records {
            writeln!(writer, "@{}", record.id)?;
            writeln!(writer, "{}", String::from_utf8_lossy(&record.sequence))?;
            writeln!(writer, "+")?;
            writeln!(writer, "{}", String::from_utf8_lossy(&record.quality))?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Check if contamination levels are within acceptable thresholds
    pub fn is_contamination_acceptable(&self, stats: &PrimerRemovalStats) -> bool {
        stats.contamination_level <= self.contamination_threshold
    }
}