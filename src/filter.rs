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
