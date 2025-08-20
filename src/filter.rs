//! Various filters to filter called variants with likely issues
use crate::model::{Call, Genotype};

/// Filter variant calls
pub trait Filter {
    /// Return true if it passes the filter
    fn filter(&self, call: &Call) -> bool;
    /// Name of the filter, used in the vcf file
    fn name(&self) -> &str;
    /// Description of the filter, used in the vcf file
    fn description(&self) -> &str;
}

// Filtering sides with too little coverage (depth)
pub struct MinCovFilter {
    description: String,
}

impl MinCovFilter {
    pub fn new(threshold: usize) -> Self {
        MinCovFilter {
            description: format!("DP < {}", threshold),
        }
    }
}

impl Filter for MinCovFilter {
    fn filter(&self, _: &Call) -> bool {
        panic!("Filtering of coverage should happen before the call")
    }

    fn name(&self) -> &str {
        "MinCov"
    }

    fn description(&self) -> &str {
        &self.description
    }
}

/// Filtering likely wrong calls via a strands bias statistic.
/// See <https://gatk.broadinstitute.org/hc/en-us/articles/360037224532-AS-StrandOddsRatio>
pub struct OddsRatioFilter {
    threshold: f64,
    description: String,
}

impl OddsRatioFilter {
    pub fn new(threshold: f64) -> Self {
        OddsRatioFilter {
            threshold,
            description: format!("Bias on strandedness. ORT >= {}", threshold),
        }
    }
}

impl Filter for OddsRatioFilter {
    fn filter(&self, call: &Call) -> bool {
        call.genotype == Genotype::HomozygousReference
            || call.genotype == Genotype::HomozygousAlternate
            || call.bias.odds_ratio < self.threshold
    }

    fn name(&self) -> &str {
        "StrandBiasOddsRatio"
    }
    fn description(&self) -> &str {
        &self.description
    }
}

/// Filter likely wrong calls if there are too many deletions
pub struct ProbabableDeletionFilter {
    threshold: f64,
    description: String,
}

impl ProbabableDeletionFilter {
    pub fn new(threshold: f64) -> Self {
        ProbabableDeletionFilter {
            threshold,
            description: format!("N deletions / (N deletions + M bases) >= {}", threshold),
        }
    }
}

impl Filter for ProbabableDeletionFilter {
    fn filter(&self, call: &Call) -> bool {
        call.genotype == Genotype::HomozygousReference
            || call.genotype == Genotype::HomozygousAlternate
            || (call.general.dels as f64 / (call.general.dels as usize + call.general.depth) as f64)
                < self.threshold
    }

    fn name(&self) -> &str {
        "ProbableDeletion"
    }

    fn description(&self) -> &str {
        &self.description
    }
}

/// Filter likely wrong calls if there are too many low quality bases (might be due to BAQ)
pub struct TooManyLowQualityFilter {
    threshold: f64,
    description: String,
}

impl TooManyLowQualityFilter {
    pub fn new(threshold: f64) -> Self {
        TooManyLowQualityFilter {
            threshold,
            description: format!(
                "N low quality reads / (N low quality reads + M bases) >= {}",
                threshold
            ),
        }
    }
}

impl Filter for TooManyLowQualityFilter {
    fn filter(&self, call: &Call) -> bool {
        call.genotype == Genotype::HomozygousReference
            || call.genotype == Genotype::HomozygousAlternate
            || (call.general.low_quals as f64
                / (call.general.low_quals as usize + call.general.depth) as f64)
                < self.threshold
    }

    fn name(&self) -> &str {
        "TooManyLowQualityReads"
    }

    fn description(&self) -> &str {
        &self.description
    }
}

pub fn failed_filters<'a>(call: &Call, filters: &'a [Box<dyn Filter>]) -> Vec<&'a dyn Filter> {
    filters
        .iter()
        .filter_map(|f| {
            if !f.filter(call) {
                Some(f.as_ref())
            } else {
                None
            }
        })
        .collect()
}
