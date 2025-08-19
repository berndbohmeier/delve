use std::collections::HashMap;

use argmin::core::{CostFunction, Error};
use argmin::core::{Executor, State};
use argmin::solver::brent::BrentOpt;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaseRead {
    pub base: u8,
    pub error: f64,
    pub strand: Strand,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Genotype {
    HomozygousReference,
    Heterozygous,
    HomozygousAlternate,
    Unknown,
}

struct BiallelicData {
    pub reference: Vec<BaseRead>,
    pub alt: Vec<BaseRead>,
}

// Filter data for the two most common alleles, sort so reference comes first
fn filter_biallelic(data: &[BaseRead], reference: u8) -> BiallelicData {
    let mut counter: HashMap<u8, usize> = HashMap::new();
    for d in data {
        *counter.entry(d.base).or_insert(0) += 1;
    }
    let mut most_common: Vec<u8> = counter.keys().copied().collect();
    most_common.sort_by_key(|b| -(counter[b] as isize));
    most_common.truncate(2);
    if !most_common.contains(&reference) {
        if most_common.len() < 2 {
            most_common.push(reference);
        } else {
            most_common[1] = reference;
        }
    }
    let mut res_reference = Vec::new();
    let mut res_alt = Vec::new();
    data.iter()
        .filter(|x| most_common.contains(&x.base))
        .for_each(|x| {
            if x.base == reference {
                res_reference.push(*x);
            } else {
                res_alt.push(*x);
            }
        });
    BiallelicData {
        reference: res_reference,
        alt: res_alt,
    }
}

// Returns (errors_matching_ref, errors_mismatching_ref)
fn error_probabilities(data: &BiallelicData) -> (Vec<f64>, Vec<f64>) {
    (
        data.reference.iter().map(|x| x.error).collect(),
        data.alt.iter().map(|x| x.error).collect(),
    )
}

// Log likelihood for a given p
fn p_log_likelihood(p: f64, errors_matching_ref: &[f64], errors_mismatching_ref: &[f64]) -> f64 {
    let mut sum = 0.0;
    for &e in errors_matching_ref {
        sum += (p * e + (1.0 - p) * (1.0 - e)).ln();
    }
    for &e in errors_mismatching_ref {
        sum += (p * (1.0 - e) + (1.0 - p) * e).ln();
    }
    sum
}

struct LogLikelihood<'a> {
    errors_matching_ref: &'a [f64],
    errors_mismatching_ref: &'a [f64],
}

impl<'a> CostFunction for LogLikelihood<'a> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        // Minimize negative log likelihood
        Ok(-p_log_likelihood(
            *p,
            self.errors_matching_ref,
            self.errors_mismatching_ref,
        ))
    }
}

// Maximize log likelihood numerically using Brent's method from argmin
fn p_max_log_likelihood(
    errors_matching_ref: &[f64],
    errors_mismatching_ref: &[f64],
    bounds: (f64, f64),
) -> f64 {
    let op = LogLikelihood {
        errors_matching_ref,
        errors_mismatching_ref,
    };
    let solver = BrentOpt::new(bounds.0, bounds.1);
    let res = Executor::new(op, solver)
        .configure(|state| state.param((bounds.0 + bounds.1) / 2.0).max_iters(100))
        .run();
    match res {
        Ok(r) => *r.state().get_best_param().unwrap(),
        Err(_) => (bounds.0 + bounds.1) / 2.0, // TODO return error
    }
}

fn log_likelihood_ratio(ll_p_ml: f64, ll_p0: f64) -> f64 {
    -2.0 * (ll_p0 - ll_p_ml)
}

struct DP4Table {
    pub forward_ref: usize,
    pub reverse_ref: usize,
    pub forward_alt: usize,
    pub reverse_alt: usize,
}

fn dp4_table(data: &BiallelicData) -> DP4Table {
    let mut table = DP4Table {
        forward_ref: 0,
        reverse_ref: 0,
        forward_alt: 0,
        reverse_alt: 0,
    };
    for x in &data.reference {
        if x.strand == Strand::Forward {
            table.forward_ref += 1;
        } else {
            table.reverse_ref += 1;
        }
    }
    for x in &data.alt {
        if x.strand == Strand::Forward {
            table.forward_alt += 1;
        } else {
            table.reverse_alt += 1;
        }
    }
    table
}

fn odds_ratio_test(dp4: &DP4Table) -> f64 {
    let r = ((dp4.forward_ref + 1) as f64 * (dp4.reverse_alt + 1) as f64)
        / ((dp4.forward_alt + 1) as f64 * (dp4.reverse_ref + 1) as f64);
    r + 1.0 / r
}

#[derive(Clone, Debug, PartialEq)]
pub struct Position {
    pub name: String,
    pub pos: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PileUpColumn {
    pub position: Position,
    pub reference: u8,
    pub data: Vec<BaseRead>,
}

pub struct Model {
    pub vaf_threshold: f64,
    pub lambda_threshold_ref: f64,
    pub lambda_thereshold_alt: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stats {
    pub mvaf: f64,
    pub vaf: f64,
    pub log_likelihood_ratio_ref: f64,
    pub log_likelihood_ratio_alt: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiasStats {
    pub odds_ratio: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct General {
    pub depth: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Call {
    pub position: Position,
    pub genotype: Genotype,
    pub ref_al: u8,
    pub alt_al: Vec<u8>,
    pub stats: Stats,
    pub bias: BiasStats,
    pub general: General,
}

impl Model {
    pub fn call(&self, col: &PileUpColumn) -> Call {
        let reference = col.reference;
        let data = &col.data;

        // Filter biallelic data
        let filtered_data = filter_biallelic(data, reference);

        let dp4 = dp4_table(&filtered_data);
        let odds_ratio = odds_ratio_test(&dp4);

        // Calculate error probabilities
        let (errors_matching_ref, errors_mismatching_ref) = error_probabilities(&filtered_data);

        let stats = self.calc(&errors_matching_ref, &errors_mismatching_ref);
        let genotype = self.genotype(&stats);

        let bias = BiasStats { odds_ratio };
        let general = General { depth: data.len() };
        Call {
            position: col.position.clone(),
            genotype,
            stats,
            bias,
            general,
            ref_al: col.reference,
            alt_al: alt_allel(&filtered_data, &genotype),
        }
    }

    fn calc(&self, errors_matching_ref: &[f64], errors_mismatching_ref: &[f64]) -> Stats {
        if errors_mismatching_ref.is_empty() {
            return Stats {
                mvaf: 0.0,
                vaf: 0.0,
                log_likelihood_ratio_ref: 0.0,
                log_likelihood_ratio_alt: f64::MAX,
            };
        }
        if errors_matching_ref.is_empty() {
            return Stats {
                mvaf: 1.0,
                vaf: 1.0,
                log_likelihood_ratio_ref: f64::MAX,
                log_likelihood_ratio_alt: 0.0,
            };
        }

        let vaf = errors_mismatching_ref.len() as f64
            / (errors_matching_ref.len() + errors_mismatching_ref.len()) as f64;
        let p_h1 = p_max_log_likelihood(errors_matching_ref, errors_mismatching_ref, (0.0, 1.0));

        if p_h1 < self.vaf_threshold {
            return Stats {
                mvaf: p_h1,
                vaf,
                log_likelihood_ratio_ref: 0.0,
                log_likelihood_ratio_alt: f64::MAX,
            };
        }

        if p_h1 > 1.0 - self.vaf_threshold {
            return Stats {
                mvaf: p_h1,
                vaf,
                log_likelihood_ratio_ref: f64::MAX,
                log_likelihood_ratio_alt: 0.0,
            };
        }

        // Compute maximum likelihood under H0, H1
        let p_ref_h0 = p_max_log_likelihood(
            errors_matching_ref,
            errors_mismatching_ref,
            (0.0, self.vaf_threshold),
        );

        let p_alt_h0 = p_max_log_likelihood(
            errors_matching_ref,
            errors_mismatching_ref,
            (1.0 - self.vaf_threshold, 1.0),
        );

        // Calculate log likelihoods
        let ll_p_h1 = p_log_likelihood(p_h1, errors_matching_ref, errors_mismatching_ref);
        let ll_p_ref_h0 = p_log_likelihood(p_ref_h0, errors_matching_ref, errors_mismatching_ref);
        let ll_p_alt_h0 = p_log_likelihood(p_alt_h0, errors_matching_ref, errors_mismatching_ref);

        // Get LRT statistics
        let lambda_ref = log_likelihood_ratio(ll_p_h1, ll_p_ref_h0);
        let lambda_alt = log_likelihood_ratio(ll_p_h1, ll_p_alt_h0);

        Stats {
            mvaf: p_h1,
            vaf,
            log_likelihood_ratio_ref: lambda_ref,
            log_likelihood_ratio_alt: lambda_alt,
        }
    }

    pub fn genotype(&self, stats: &Stats) -> Genotype {
        if stats.log_likelihood_ratio_ref < self.lambda_threshold_ref {
            return Genotype::HomozygousReference;
        }
        if stats.log_likelihood_ratio_alt < self.lambda_thereshold_alt {
            return Genotype::HomozygousAlternate;
        }
        Genotype::Heterozygous
    }
}

// Return the alt allel to report
fn alt_allel(biallelic: &BiallelicData, genotype: &Genotype) -> Vec<u8> {
    match genotype {
        Genotype::HomozygousReference | Genotype::Unknown => Vec::new(),
        Genotype::HomozygousAlternate | Genotype::Heterozygous => vec![
            biallelic
                .alt
                .last()
                .expect("should have an alt for alt calls")
                .base,
        ],
    }
}

// write tests for the model
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_filter_biallelic() {
        let data = vec![
            BaseRead {
                base: b'A',
                error: 0.01,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'A',
                error: 0.02,
                strand: Strand::Reverse,
            },
            BaseRead {
                base: b'C',
                error: 0.03,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'T',
                error: 0.04,
                strand: Strand::Reverse,
            },
        ];
        let reference = b'T';
        let biallelic_data = filter_biallelic(&data, reference);
        assert_eq!(biallelic_data.reference.len(), 1);
        assert_eq!(biallelic_data.alt.len(), 2);
    }
    #[test]
    fn test_filter_biallelic_no_reference() {
        let data = vec![
            BaseRead {
                base: b'A',
                error: 0.01,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'A',
                error: 0.02,
                strand: Strand::Reverse,
            },
            BaseRead {
                base: b'C',
                error: 0.03,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'T',
                error: 0.04,
                strand: Strand::Reverse,
            },
        ];
        let reference = b'G';
        let filtered = filter_biallelic(&data, reference);
        assert_eq!(filtered.reference.len(), 0);
        assert_eq!(filtered.alt.len(), 2);
    }

    #[test]
    fn test_dp4_table() {
        let data = vec![
            BaseRead {
                base: b'A',
                error: 0.01,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'A',
                error: 0.02,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'T',
                error: 0.03,
                strand: Strand::Forward,
            },
            BaseRead {
                base: b'T',
                error: 0.04,
                strand: Strand::Reverse,
            },
        ];
        let reference = b'T';
        let table = dp4_table(&filter_biallelic(&data, reference));
        assert_eq!(table.forward_ref, 1);
        assert_eq!(table.reverse_ref, 1);
        assert_eq!(table.forward_alt, 2);
        assert_eq!(table.reverse_alt, 0);
    }
}
