//! Variant caller for mixed infections, mainly with the goal to handle Malaria infections
use std::iter::once;

use anyhow::Ok;
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use rust_htslib::bam;
use rust_htslib::htslib;

use crate::filter::Filter;
use crate::io::Region;
use crate::io::{VCFWriter, load_regions, parse_region};

mod filter;
mod io;
mod model;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Call variants from a BAM file
    Call(CallArgs),
}

/// Arguments for variant calling
#[derive(Parser, Debug)]
struct CallArgs {
    /// BAM file
    bamfile: String,

    /// Reference FASTA file
    #[clap(
        short = 'f',
        long = "fasta-ref",
        value_name = "FASTA",
        required = true,
        value_name = "FILE"
    )]
    fasta_ref: String,

    /// Regions file
    #[clap(short = 'R', long = "regions-file", value_name = "FILE")]
    regions_file: Option<String>,

    /// Region string
    #[clap(short = 'r', long = "region")]
    region: Option<String>,

    /// Sample name
    #[clap(
        short = 's',
        long = "sample_name",
        default_value = "sample",
        value_name = "NAME"
    )]
    sample_name: String,

    /// Output file
    #[clap(short = 'o', long = "output", value_name = "FILE")]
    output: Option<String>,

    /// Minimum mapping quality
    #[clap(short = 'q', long = "min-MQ", default_value_t = 0)]
    min_mq: i32,

    /// Minimum base quality
    #[clap(short = 'Q', long = "min-BQ", default_value_t = 20)]
    min_bq: i32,

    /// Minimum coverage
    #[clap(long = "min-cov", default_value_t = 10)]
    min_cov: i32,

    /// Maximum coverage
    #[clap(long = "max-cov", default_value_t = 5000)]
    max_cov: i32,

    /// Minimum VAF
    // #[clap(long = "min-vaf", default_value_t = 0.01)]
    // min_vaf: f64,

    /// Truncate regions
    #[clap(long = "truncate-regions", default_value_t = 0, value_name = "INT")]
    truncate_regions: usize,

    /// Compute BAQ
    #[clap(long = "compute-baq")]
    compute_baq: bool,

    /// Strand bias odds ratio
    #[clap(
        long = "strand-bias-odds-ratio",
        default_value_t = 7.0,
        value_name = "FLOAT"
    )]
    strand_bias_odds_ratio: f64,

    /// Deletion filter threshold. Positions with higher ratio of deletions will be filtered.
    #[clap(
        long = "deletion-filter-threshold",
        default_value_t = 0.8,
        value_name = "FLOAT"
    )]
    deletion_filter_threshold: f64,

    /// Too many low quality reads filter threshold. Positions with higher ratio of low quality reads will be filtered.
    #[clap(
        long = "low-qual-reads-filter-threshold",
        default_value_t = 0.8,
        value_name = "FLOAT"
    )]
    low_quality_reads_filter_threshold: f64,

    /// Model parameters (comma-separated floats)
    #[clap(
        long = "model-params",
        use_value_delimiter = true,
        value_delimiter = ',',
        default_value = "8.0,8.0,0.01",
        value_name = "LRT_REF,LRT_ALT,H0_VAF"
    )]
    model_params: Vec<f64>,

    /// Show only variants
    #[clap(short = 'v', long = "variants-only")]
    variants_only: bool,

    /// Apply filters
    #[clap(long = "apply-filters", value_name = "LIST")]
    apply_filters: Option<String>,

    /// Set genotypes of failed samples to missing value (.) or reference (0)
    #[clap(long = "set-failed-GTs", value_name = "TYPE", value_parser = ["0", "."])]
    set_failed_gts: Option<String>,
}

fn call_variants(args: CallArgs) -> Result<()> {
    let fasta_reader = rust_htslib::faidx::Reader::from_path(&args.fasta_ref)
        .with_context(|| format!("failed to open FASTA file `{}`", &args.fasta_ref))?;

    let seq_names = fasta_reader
        .seq_names()
        .with_context(|| "failed to get sequence names from FASTA file")?;
    let contigs: Vec<_> = seq_names
        .iter()
        .map(|name| (name.as_str(), fasta_reader.fetch_seq_len(name) as usize))
        .collect();

    let regions = if let Some(region) = args.region {
        vec![parse_region(&region).with_context(|| format!("failed to parse region `{region}`"))?]
    } else if let Some(regions_file) = args.regions_file {
        load_regions(&regions_file)
            .with_context(|| format!("failed to load regions from file `{regions_file}`"))?
    } else {
        contigs
            .iter()
            .map(|(name, len)| Region {
                chrom: name.to_string(),
                start: 0,
                end: *len,
            })
            .collect()
    };

    let filters: Vec<Box<dyn filter::Filter>> = vec![
        Box::new(filter::OddsRatioFilter::new(args.strand_bias_odds_ratio)),
        Box::new(filter::TooManyLowQualityFilter::new(
            args.low_quality_reads_filter_threshold,
        )),
        Box::new(filter::ProbabableDeletionFilter::new(
            args.deletion_filter_threshold,
        )),
    ];

    let min_cov_filter = filter::MinCovFilter::new(args.min_cov as usize);

    let all_filters = filters
        .iter()
        .map(|f| f.as_ref())
        .chain(once(&min_cov_filter as &dyn Filter))
        .collect::<Vec<_>>();
    let mut vcf_writer = VCFWriter::create(
        args.output.as_ref().map_or("-", |o| o.as_str()),
        &contigs,
        &args.sample_name,
        &all_filters,
        true,
    )
    .with_context(|| "failed to create vcf file")?;
    let mut bam_reader = bam::IndexedReader::from_path(&args.bamfile)
        .with_context(|| format!("failed to open BAM file `{}`", &args.bamfile))?;
    bam_reader.set_filters(
        htslib::BAM_FUNMAP
            | htslib::BAM_FSECONDARY
            | htslib::BAM_FQCFAIL
            | htslib::BAM_FSUPPLEMENTARY
            | htslib::BAM_FDUP,
    );

    if args.compute_baq {
        bam_reader
            .enable_baq_computation(
                &args.fasta_ref,
                (htslib::htsRealnFlags_BAQ_APPLY
                    | htslib::htsRealnFlags_BAQ_EXTEND
                    | htslib::htsRealnFlags_BAQ_ONT) as usize,
            )
            .with_context(|| "failed to enable BAQ computation")?;
    }
    for region in regions {
        let options = io::PileUpOptions {
            min_mq: args.min_mq as u8,
            min_bq: args.min_bq as u8,
            max_cov: args.max_cov as u32,
            truncate_regions: args.truncate_regions,
        };
        let iter = io::pileup_region(&mut bam_reader, &fasta_reader, &region, &options)
            .with_context(|| "failed to create pileup iterator")?;

        let model = model::Model {
            vaf_threshold: args.model_params[2],
            lambda_threshold_ref: args.model_params[0],
            lambda_thereshold_alt: args.model_params[1],
        };

        let apply_filters: Vec<_> = args
            .apply_filters
            .as_ref()
            .map_or_else(Vec::new, |s| s.split(",").collect());

        for column in iter {
            let column = column.with_context(|| "failed to pileup column")?;
            let (mut call, failed_filters) = if column.data.len() < args.min_cov as usize {
                let call = model::Call {
                    position: column.position.clone(),
                    genotype: model::Genotype::Unknown,
                    bias: model::BiasStats { odds_ratio: 0.0 },
                    ref_al: column.reference,
                    alt_al: vec![],
                    stats: model::Stats {
                        mvaf: 0.,
                        vaf: 0.,
                        log_likelihood_ratio_ref: 0.,
                        log_likelihood_ratio_alt: 0.,
                    },
                    general: model::General {
                        dels: column.dels,
                        low_quals: column.low_quals,
                        depth: column.data.len(),
                    },
                };
                let failed_filters: Vec<&dyn Filter> = vec![&min_cov_filter];
                (call, failed_filters)
            } else {
                let call = model.call(&column);
                let failed_filters = filter::failed_filters(&call, &filters);
                (call, failed_filters)
            };
            if args.variants_only
                && (call.genotype != model::Genotype::Heterozygous
                    && call.genotype != model::Genotype::HomozygousAlternate)
            {
                // Filter out non variant calls
                continue;
            }
            if io::should_filter(
                &apply_filters,
                &failed_filters.iter().map(|f| f.name()).collect::<Vec<_>>(),
            ) {
                // Filter out sites with no matching filter
                continue;
            }
            if args.set_failed_gts.is_some() && !failed_filters.is_empty() {
                match args.set_failed_gts.as_deref() {
                    Some(".") => {
                        call.genotype = model::Genotype::Unknown;
                    }
                    Some("0") => {
                        if !failed_filters.is_empty() {
                            call.genotype = model::Genotype::HomozygousReference;
                        }
                    }
                    Some(_) => {
                        return Err(anyhow::anyhow!(
                            "unknown value for --set-failed-GTs, use . or 0"
                        ));
                    }
                    None => {}
                }
            }
            vcf_writer
                .write_call(&call, &failed_filters)
                .with_context(|| "failed to write call to VCF")?;
        }
    }
    Ok(())
}

fn main() {
    let cli = Cli::parse();

    if let Err(error) = match cli.command {
        Commands::Call(args) => call_variants(args),
    } {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert();
}
