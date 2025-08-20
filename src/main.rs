use clap::Parser;
use rust_htslib::bam;
use rust_htslib::htslib;
use std::iter::once;

use crate::io::{VCFWriter, load_regions, parse_region};

mod filter;
mod io;
mod model;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
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
}

fn main() {
    let cli = Cli::parse();

    let regions = if let Some(region) = cli.region {
        vec![parse_region(&region).expect("Failed to parse region")]
    } else if let Some(regions_file) = cli.regions_file {
        load_regions(&regions_file).expect("Failed to load regions")
    } else {
        panic!("Needs one of region or region file");
    };

    let fasta_reader =
        rust_htslib::faidx::Reader::from_path(&cli.fasta_ref).expect("Failed to open FASTA file");

    let seq_names = fasta_reader
        .seq_names()
        .expect("Failed to get sequence names from FASTA file");
    let contigs: Vec<_> = seq_names
        .iter()
        .map(|name| (name.as_str(), fasta_reader.fetch_seq_len(name) as usize))
        .collect();

    let filters: Vec<Box<dyn filter::Filter>> = vec![
        Box::new(filter::OddsRatioFilter::new(cli.strand_bias_odds_ratio)),
        Box::new(filter::TooManyLowQualityFilter::new(
            cli.low_quality_reads_filter_threshold,
        )),
        Box::new(filter::ProbabableDeletionFilter::new(
            cli.deletion_filter_threshold,
        )),
    ];

    let filter_names = filters
        .iter()
        .map(|f| f.name())
        .chain(once("MinCov"))
        .collect::<Vec<_>>();
    let mut vcf_writer = VCFWriter::new(
        cli.output.as_ref().map_or("-", |o| o.as_str()),
        &contigs,
        &cli.sample_name,
        &filter_names,
        true,
    )
    .expect("Failed to create VCF writer");
    let mut bam_reader =
        bam::IndexedReader::from_path(&cli.bamfile).expect("Failed to open BAM file");
    bam_reader.set_filters(
        htslib::BAM_FUNMAP
            | htslib::BAM_FSECONDARY
            | htslib::BAM_FQCFAIL
            | htslib::BAM_FSUPPLEMENTARY
            | htslib::BAM_FDUP,
    );

    if cli.compute_baq {
        bam_reader
            .enable_baq_computation(
                &cli.fasta_ref,
                (htslib::htsRealnFlags_BAQ_APPLY
                    | htslib::htsRealnFlags_BAQ_EXTEND
                    | htslib::htsRealnFlags_BAQ_ONT) as usize,
            )
            .expect("Faled to enable baq computation");
    }
    for region in regions {
        let options = io::PileUpOptions {
            min_mq: cli.min_mq as u8,
            min_bq: cli.min_bq as u8,
            max_cov: cli.max_cov as u32,
            truncate_regions: cli.truncate_regions,
        };
        let iter = io::pileup_region(&mut bam_reader, &fasta_reader, &region, &options)
            .expect("Failed to create pileup iterator");

        let model = model::Model {
            vaf_threshold: cli.model_params[2],
            lambda_threshold_ref: cli.model_params[0],
            lambda_thereshold_alt: cli.model_params[1],
        };

        let apply_filters: Vec<_> = cli
            .apply_filters
            .as_ref()
            .map_or_else(Vec::new, |s| s.split(",").collect());

        for column in iter {
            let (call, failed_filters) = if column.data.len() < cli.min_cov as usize {
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
                let failed_filters = vec!["MinCov"];
                (call, failed_filters)
            } else {
                let call = model.call(&column);
                let failed_filters = filter::failed_filters(&call, &filters);
                (call, failed_filters)
            };
            if cli.variants_only
                && (call.genotype != model::Genotype::Heterozygous
                    && call.genotype != model::Genotype::HomozygousAlternate)
            {
                // Filter out non variant calls
                continue;
            }
            if io::should_filter(&apply_filters, &failed_filters) {
                // Filter out sites with no matching filter
                continue;
            }
            vcf_writer
                .write_call(&call, &failed_filters)
                .expect("Failed to write call to VCF");
        }
    }
}

#[test]
fn verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert();
}
