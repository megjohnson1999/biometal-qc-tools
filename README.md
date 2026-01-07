# Biometal QC Tools

Fast QC tools for virome analysis using proven biometal primitives.

## Overview

This project provides high-performance quality control tools for viral metagenomics (virome) analysis, built on top of the proven biometal primitives library. These tools are designed to replace slower external QC tools with fast, composition-based alternatives optimized for virome workflows.

## Tools Included

### 1. biometal-quality-stats (Priority)
Fast quality assessment tool to replace FastQC functionality:
- Base composition analysis using `base_counting`
- GC content calculation using `gc_content`
- Quality score distributions using `quality_filter`
- Sequence complexity assessment using `complexity`

### 2. biometal-contamination-screen
PhiX and vector contamination detection:
- Pattern matching for known contamination sequences
- Composition-based contamination identification

### 3. biometal-vlp-assessment
VLP (Virus-Like Particle) success metrics:
- Composition-based assessment (not alignment-based like ViromeQC)
- GC distribution scoring
- Sequence diversity and complexity analysis
- Compositional evenness metrics

### 4. biometal-qc-summary
Multi-sample QC reporting:
- Aggregate statistics across samples
- Pass/fail determination based on configurable thresholds
- JSON output for downstream analysis

## Proven Biometal Primitives Used

This project uses only the **proven, stable biometal primitives**:
- ✅ `base_counting` - Base composition analysis
- ✅ `gc_content` - GC content calculation
- ✅ `quality_filter` - Quality score filtering
- ✅ `complexity` - Sequence complexity assessment
- ✅ `pattern_match` - Pattern matching for contamination
- ✅ `trimming` - Sequence trimming operations

**Note**: Alignment primitives are intentionally excluded as they still have quality issues.

## Installation

```bash
# Clone the repository
git clone git@github.com:megjohnson1999/biometal-qc-tools.git
cd biometal-qc-tools

# Build all tools
cargo build --release

# Tools will be available in target/release/
```

## Usage

### Quality Statistics
```bash
./target/release/biometal-quality-stats -i sample.fastq -o quality_stats.json
```

### Contamination Screening
```bash
./target/release/biometal-contamination-screen -i sample.fastq -o contamination_report.json
```

### VLP Assessment
```bash
./target/release/biometal-vlp-assessment -i sample.fastq -o vlp_assessment.json
```

### Multi-Sample Summary
```bash
./target/release/biometal-qc-summary -i qc_results_dir/ -o qc_summary.json
```

## Dependencies

- **biometal**: Local dependency from `../biometal` (proven primitives only)
- **clap**: Command-line argument parsing
- **serde**: JSON serialization/deserialization
- **anyhow**: Error handling

## Project Goals

1. **Performance**: Replace slow external tools with fast biometal-based alternatives
2. **Reliability**: Use only proven, stable biometal primitives
3. **Practicality**: Build tools that directly solve real virome QC needs
4. **Integration**: Easy integration into existing virome analysis pipelines

## License

MIT OR Apache-2.0