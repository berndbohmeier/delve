//! IO features to read and write various file types
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Error};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Writer};
use rust_htslib::htslib::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP};
use rust_htslib::{bam, bam::Read, faidx};
use thiserror::Error;

use crate::filter::Filter;
use crate::model::{BaseRead, Call, Genotype, PileUpColumn, Position, Strand};

// struct to represent a region in the genome
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
}

pub struct PileUpOptions {
    pub min_mq: u8,
    pub min_bq: u8,
    pub max_cov: u32,
}

const DEFAULT_FLAG_FILTER: u32 =
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;

pub fn pileup_region(
    bam_reader: &mut bam::IndexedReader,
    fasta_reader: &faidx::Reader,
    region: &Region,
    options: &PileUpOptions,
) -> Result<impl Iterator<Item = Result<PileUpColumn, Error>>, Error> {
    bam_reader
        .fetch((
            region.chrom.as_str(),
            region.start as i64,
            region.end as i64,
        ))
        .with_context(|| "failed to fetch region")?;

    // Get the reference base for the regio
    let reference_region = fasta_reader
        .fetch_seq(region.chrom.as_str(), region.start, region.end)
        .with_context(|| "failed to fetch reference region")?;

    let mut pileup_iter = bam_reader.pileup();
    pileup_iter.set_max_depth(options.max_cov);
    Ok(pileup_iter.filter_map(move |pileup| {
        let p = match pileup {
            Ok(p) => p,
            Err(e) => return Some(Err(Error::new(e))),
        };
        if (p.pos() as usize) < region.start || (p.pos() as usize) >= region.end {
            // Ignore reads outside of the region
            return None;
        }
        let mut dels = 0;
        let mut low_quals = 0;
        let data: Vec<BaseRead> = p
            .alignments()
            .filter_map(|alignment| {
                if alignment.is_del() {
                    dels += 1;
                    return None;
                }
                if alignment.is_refskip() {
                    return None;
                }
                let record = alignment.record();
                if record.flags() & DEFAULT_FLAG_FILTER as u16 != 0 {
                    return None; // Skip filtered reads
                }
                if options.min_mq > 0 && record.mapq() < options.min_mq {
                    low_quals += 1;
                    return None; // Skip low mapping quality reads
                }
                assert_eq!(record.seq_len(), record.qual().len());
                let base = record.seq()[alignment.qpos().unwrap()];
                let qual = record.qual()[alignment.qpos().unwrap()];
                if qual < options.min_bq {
                    low_quals += 1;
                    return None; // Skip low quality bases
                }
                let strand = if record.is_reverse() {
                    Strand::Reverse
                } else {
                    Strand::Forward
                };
                Some(BaseRead {
                    base,
                    error: 10f64.powf(-(qual as f64) / 10.0),
                    strand,
                })
            })
            .collect();

        Some(Ok(PileUpColumn {
            position: Position {
                name: region.chrom.clone(),
                pos: p.pos() as usize,
            },
            reference: reference_region[p.pos() as usize - region.start],
            data,
            dels,
            low_quals,
        }))
    }))
}

#[derive(Error, Debug)]
pub enum LoadRegionsError<'a> {
    #[error("failed to read file")]
    IO(#[from] std::io::Error),
    #[error("invalid {1} at line {2}")]
    InvalidInteger(#[source] std::num::ParseIntError, &'a str, usize),
    #[error("missing {0} in line {1}")]
    MissingField(&'a str, usize),
}

// Loads regions from a file: chrom, start, end (tab or space separated)
pub fn load_regions(filename: &str) -> Result<Vec<Region>, LoadRegionsError<'static>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut result = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let mut parts = line.split_whitespace();
        let chrom = parts
            .next()
            .ok_or(LoadRegionsError::MissingField("chrom", i + 1))?
            .to_string();
        let start = parts
            .next()
            .ok_or(LoadRegionsError::MissingField("start", i + 1))?
            .parse::<usize>()
            .map_err(|err| LoadRegionsError::InvalidInteger(err, "start", i + 1))?;
        let end = parts
            .next()
            .ok_or(LoadRegionsError::MissingField("end", i + 1))?
            .parse::<usize>()
            .map_err(|err| LoadRegionsError::InvalidInteger(err, "end", i + 1))?;
        result.push(Region { chrom, start, end });
    }
    Ok(result)
}

#[derive(Error, Debug)]
pub enum ParseRegionError<'a> {
    #[error("invalid integer")]
    InvalidInteger(#[from] std::num::ParseIntError),
    #[error("missing {0}")]
    MissingField(&'a str),
}

// Parses a region string like "chr1:100-200" into (chrom, start, end)
pub fn parse_region(region_str: &str) -> Result<Region, ParseRegionError<'static>> {
    let mut parts = region_str.split(':');
    let chrom = parts
        .next()
        .ok_or(ParseRegionError::MissingField("chrom"))?
        .to_string();
    let range = parts
        .next()
        .ok_or(ParseRegionError::MissingField("range"))?;
    let mut range_parts = range.split('-');
    let start = range_parts
        .next()
        .ok_or(ParseRegionError::MissingField("start"))?
        .parse::<usize>()?;
    let end = range_parts
        .next()
        .ok_or(ParseRegionError::MissingField("end"))?
        .parse::<usize>()?;
    Ok(Region { chrom, start, end })
}

fn format_header(id: &str, number: usize, type_: &str, description: &str) -> Vec<u8> {
    format!(
        r#"##FORMAT=<ID={id},Number={number},Type={type_},Description="{description}">"#,
        id = id,
        number = number,
        type_ = type_,
        description = description
    )
    .into_bytes()
}

#[derive(Error, Debug)]
#[error("vcf write error")]
pub struct VCFError(#[from] rust_htslib::errors::Error);

pub struct VCFWriter {
    vcf: Writer,
}

impl VCFWriter {
    pub fn create(
        path: &str,
        contigs: &[(&str, usize)],
        sample_name: &str,
        filters: &[&dyn Filter],
        compressed: bool,
    ) -> Result<Self, VCFError> {
        let mut header = Header::new();

        header.push_record(format!("##delveVersion={}", env!("CARGO_PKG_VERSION")).as_bytes());
        header.push_record(
            format!(
                "##delveCommand={}",
                env::args().skip(1).collect::<Vec<_>>().join(" ")
            )
            .as_bytes(),
        );

        for (name, length) in contigs {
            header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
        }

        for filter in filters {
            header.push_record(
                format!(
                    "##FILTER=<ID={},Description=\"{}\">",
                    filter.name(),
                    filter.description()
                )
                .as_bytes(),
            );
        }

        header.push_record(&format_header("GT", 1, "String", "Genotype"));
        header.push_record(&format_header("DP", 1, "Integer", "Read depth"));
        header.push_record(&format_header("ORT", 1, "Float", "Odds ratio test value"));
        header.push_record(&format_header(
            "VAF",
            1,
            "Float",
            "Variant allele frequency",
        ));
        header.push_record(&format_header("MVAF", 1, "Float", "Model estimated VAF"));
        header.push_record(&format_header(
            "LRT",
            2,
            "Float",
            "LRT statistic againts p=0 and p=1",
        ));
        header.push_record(&format_header("DEL", 1, "Integer", "Number of deletions"));
        header.push_record(&format_header(
            "LQR",
            1,
            "Integer",
            "Number of low quality reads",
        ));

        header.push_sample(sample_name.as_bytes());
        let vcf = Writer::from_path(path, &header, compressed, Format::Vcf)?;
        Ok(VCFWriter { vcf })
    }

    pub fn write_call(
        &mut self,
        call: &Call,
        failed_filters: &[&dyn Filter],
    ) -> Result<(), VCFError> {
        let mut record = self.vcf.empty_record();
        let header = self.vcf.header();
        record.set_rid(header.name2rid(call.position.name.as_bytes()).ok());
        record.set_pos(call.position.pos as i64);
        record.set_alleles(&[
            &[call.ref_al],
            if !call.alt_al.is_empty() {
                &call.alt_al
            } else {
                b"."
            },
        ])?;
        let alleles = match call.genotype {
            Genotype::HomozygousReference => {
                &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
            }
            Genotype::Heterozygous => &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)],
            Genotype::HomozygousAlternate => {
                &[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]
            }
            Genotype::Unknown => &[
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ],
        };
        record.push_genotypes(alleles)?;
        record.push_format_integer(b"DP", &[call.general.depth as i32])?;
        record.push_format_integer(b"DEL", &[call.general.dels as i32])?;
        record.push_format_integer(b"LQR", &[call.general.low_quals as i32])?;
        record.push_format_float(b"ORT", &[call.bias.odds_ratio as f32])?;
        record.push_format_float(b"VAF", &[call.stats.vaf as f32])?;
        record.push_format_float(b"MVAF", &[call.stats.mvaf as f32])?;
        record.push_format_float(
            b"LRT",
            &[
                call.stats.log_likelihood_ratio_ref as f32,
                call.stats.log_likelihood_ratio_alt as f32,
            ],
        )?;
        let filter_ids: Vec<_> = if !failed_filters.is_empty() {
            failed_filters
                .iter()
                .map(|f| {
                    header
                        .name_to_id(f.name().as_bytes())
                        .expect("should have filter {} in header")
                })
                .collect()
        } else {
            vec![
                header
                    .name_to_id("PASS".as_bytes())
                    .expect("should have a pass filter"),
            ]
        };
        record.set_filters(&filter_ids.iter().collect::<Vec<_>>())?;
        self.vcf.write(&record)?;
        Ok(())
    }
}

pub fn should_filter(apply_filters: &[&str], failed_filters: &[&str]) -> bool {
    !apply_filters.is_empty()
        && !apply_filters.iter().any(|filter| match *filter {
            "PASS" => failed_filters.is_empty(),
            other => failed_filters.contains(&other),
        })
}
