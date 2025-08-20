use crate::model::{Call, Genotype};

pub trait Filter {
    fn filter(&self, call: &Call) -> bool;

    fn name(&self) -> &str;
}

pub struct OddsRatioFilter {
    threshold: f64,
}

impl OddsRatioFilter {
    pub fn new(threshold: f64) -> Self {
        OddsRatioFilter { threshold }
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
}
pub struct ProbabableDeletionFilter {
    threshold: f64,
}

impl ProbabableDeletionFilter {
    pub fn new(threshold: f64) -> Self {
        ProbabableDeletionFilter { threshold }
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
}
pub struct TooManyLowQualityFilter {
    threshold: f64,
}

impl TooManyLowQualityFilter {
    pub fn new(threshold: f64) -> Self {
        TooManyLowQualityFilter { threshold }
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
}

pub fn failed_filters<'a>(call: &Call, filters: &'a [Box<dyn Filter>]) -> Vec<&'a str> {
    filters
        .iter()
        .filter_map(|f| {
            if !f.filter(call) {
                Some(f.name())
            } else {
                None
            }
        })
        .collect()
}
