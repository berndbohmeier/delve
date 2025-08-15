use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Writer};
use rust_htslib::htslib::{BAM_FDUP, BAM_FQCFAIL, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP};
use rust_htslib::{bam, bam::Read, faidx};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

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
    pub truncate_regions: usize,
}

const DEFAULT_FLAG_FILTER: u32 =
    BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP;

pub fn pileup_region(
    bam_reader: &mut bam::IndexedReader,
    fasta_reader: &faidx::Reader,
    region: &Region,
    options: &PileUpOptions,
) -> Result<impl Iterator<Item = PileUpColumn>, Box<dyn std::error::Error>> {
    bam_reader.fetch((
        region.chrom.as_str(),
        region.start as i64,
        region.end as i64,
    ))?;

    // Get the reference base for the regio
    let reference_region = fasta_reader
        .fetch_seq(region.chrom.as_str(), region.start, region.end)
        .expect("Failed to fetch reference base"); // TODO: handle error

    let mut pileup_iter = bam_reader.pileup();
    pileup_iter.set_max_depth(options.max_cov);
    Ok(pileup_iter.filter_map(move |pileup| {
        let p = pileup.unwrap(); // TODO: handle error
        if (p.pos() as usize) < region.start + options.truncate_regions
            || (p.pos() as usize) > region.end - options.truncate_regions
        {
            // Ignore reads outside of the region
            return None;
        }
        let data: Vec<BaseRead> = p
            .alignments()
            .filter_map(|alignment| {
                if alignment.is_del() || alignment.is_refskip() {
                    return None;
                }
                let record = alignment.record();
                if record.flags() & DEFAULT_FLAG_FILTER as u16 != 0 {
                    return None; // Skip filtered reads
                }
                if options.min_mq > 0 && record.mapq() < options.min_mq {
                    return None; // Skip low mapping quality reads
                }
                assert_eq!(record.seq_len(), record.qual().len());
                if record.seq_len() == 0 {
                    // TODO, why is this needed?
                    return None; // Skip empty reads
                }
                let base = record.seq()[alignment.qpos().unwrap()];
                let qual = record.qual()[alignment.qpos().unwrap()];
                if qual < options.min_bq {
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

        Some(PileUpColumn {
            position: Position {
                name: region.chrom.clone(),
                pos: p.pos() as usize,
            },
            reference: reference_region[p.pos() as usize - region.start],
            data,
        })
    }))
}

// Loads regions from a file: chrom, start, end (tab or space separated)
pub fn load_regions(filename: &str) -> Result<Vec<Region>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut result = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split_whitespace();
        let chrom = parts.next().unwrap_or("").to_string(); // TODO error handling
        let start = parts.next().unwrap_or("0").parse::<usize>().unwrap_or(0);
        let end = parts.next().unwrap_or("0").parse::<usize>().unwrap_or(0);
        result.push(Region { chrom, start, end });
    }
    Ok(result)
}

// Parses a region string like "chr1:100-200" into (chrom, start, end)
pub fn parse_region(region_str: &str) -> Result<Region, Box<dyn std::error::Error>> {
    let mut parts = region_str.split(':');
    let chrom = parts.next().ok_or("Missing chrom")?.to_string();
    let range = parts.next().ok_or("Missing range")?;
    let mut range_parts = range.split('-');
    let start = range_parts
        .next()
        .ok_or("Missing start")?
        .parse::<usize>()?;
    let end = range_parts.next().ok_or("Missing end")?.parse::<usize>()?;
    Ok(Region { chrom, start, end })
}

pub struct VCFWriter {
    vcf: Writer,
}

impl VCFWriter {
    pub fn new(
        contigs: &[(&str, usize)],
        sample_name: &str,
        filters: &[&str],
    ) -> Result<Self, Box<dyn std::error::Error>> {
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
                format!("##FILTER=<ID={},Description=\"{}\">", filter, filter).as_bytes(),
            );
        }

        header.push_record(
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#.as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#.as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=ORT,Number=1,Type=Float,Description="Odds ratio test value">"#
                .as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">"#
                .as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=MVAF,Number=1,Type=Float,Description="Model estimated VAF">"#
                .as_bytes(),
        );
        header.push_record(
            r#"##FORMAT=<ID=LRT,Number=2,Type=Float,Description="LRT statistic againts p=0 and p=1">"#.as_bytes(),
        );

        header.push_sample(sample_name.as_bytes());
        let vcf = Writer::from_stdout(&header, true, Format::Vcf)?;
        Ok(VCFWriter { vcf })
    }

    pub fn write_call(
        &mut self,
        call: &Call,
        failed_filters: &[&str],
    ) -> Result<(), Box<dyn std::error::Error>> {
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
                        .name_to_id(f.as_bytes())
                        .expect("Filter not found in header")
                })
                .collect()
        } else {
            vec![
                header
                    .name_to_id("PASS".as_bytes())
                    .expect("Pass filter missing"),
            ]
        };
        record.set_filters(&filter_ids.iter().collect::<Vec<_>>())?;
        self.vcf.write(&record)?;
        Ok(())
    }
}
