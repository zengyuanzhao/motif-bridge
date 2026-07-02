//! homer2meme - Convert HOMER motif format to MEME motif format
//!
//! Build (Cargo, recommended - includes .gz support):
//!   cd rust_scripts && cargo build --release
//!   # Binaries: rust_scripts/target/release/meme2homer  homer2meme
//!
//! Install to PATH:
//!   cd rust_scripts && cargo install --path .
//!
//! Usage:
//!   homer2meme -i motifs.homer > motifs.meme
//!   homer2meme -i motifs.homer.gz > motifs.meme
//!   homer2meme -i motifs.json -f json > motifs.meme
//!   homer2meme -i motifs.homer --input-format logodds > motifs.meme
//!   cat motifs.homer | homer2meme -i -

use clap::{ArgAction, Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
use motif_bridge::io::{
    read_homer_with_strict, read_json_with_strict, write_meme_with_background, MatrixType,
};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};

#[derive(Parser)]
#[command(
    name = "homer2meme",
    about = "Convert HOMER motif format to MEME format.",
    version,
    after_help = "Examples:\n  homer2meme -i motifs.homer > motifs.meme\n  homer2meme -i motifs.homer.gz > motifs.meme\n  homer2meme -i motifs.homer -e \"CTCF/Jaspar\"\n  homer2meme -i motifs.json -f json > motifs.meme\n  homer2meme -i motifs.homer --input-format logodds\n  cat motifs.homer | homer2meme -i -\n"
)]
struct Args {
    /// Input HOMER motif file (.homer, .motif, .json, or .gz), or '-' for stdin
    #[arg(short = 'i', value_name = "<file>")]
    input: String,
    /// Extract only specified motif by id or description
    #[arg(short = 'e', default_value = "")]
    extract: String,
    /// Pseudocount for log-odds -> probability conversion
    #[arg(short = 'a', value_parser = parse_positive, default_value = "0.01")]
    pseudocount: f64,
    /// Background probability scalar or comma-separated vector
    #[arg(short = 'b', long = "background", default_value = "0.25")]
    background: String,
    /// Input format: homer (default) or json
    #[arg(short = 'f', long = "format", value_enum, default_value_t = InputFormat::Homer)]
    input_format: InputFormat,
    /// Matrix type: auto (default), logodds, or probability
    #[arg(long = "input-format", value_enum, default_value_t = MatrixTypeArg::Auto)]
    matrix_type: MatrixTypeArg,
    /// Alphabet for HOMER input
    #[arg(long = "alphabet", value_enum, default_value_t = Alphabet::Acgt)]
    alphabet: Alphabet,
    /// Output the reverse complement of the motif (DNA/RNA only)
    #[arg(long = "rc", action = ArgAction::SetTrue)]
    rc: bool,
    /// Trim edges with information content below threshold
    #[arg(long = "trim-edges", value_parser = parse_nonneg, default_value = "0.0")]
    trim_edges: f64,
    /// Filter out motifs with total information content below threshold
    #[arg(long = "min-ic", value_parser = parse_nonneg, default_value = "0.0")]
    min_ic: f64,
    /// Override MEME nsites metadata in output
    #[arg(long = "nsites", value_parser = parse_positive_usize)]
    nsites: Option<usize>,
    /// Override MEME E metadata in output
    #[arg(long = "evalue", value_parser = parse_nonneg)]
    evalue: Option<f64>,
    /// Renormalize each row before writing MEME output
    #[arg(long = "renormalize", action = ArgAction::SetTrue)]
    renormalize: bool,
    /// Fail on malformed matrix rows, ambiguous auto-detection, or invalid probability rows
    #[arg(long = "strict", action = ArgAction::SetTrue)]
    strict: bool,
}

#[derive(Clone, Copy, ValueEnum)]
enum InputFormat {
    #[value(name = "homer")]
    Homer,
    #[value(name = "json")]
    Json,
}

#[derive(Clone, Copy, ValueEnum)]
enum MatrixTypeArg {
    #[value(name = "auto")]
    Auto,
    #[value(name = "logodds")]
    Logodds,
    #[value(name = "probability")]
    Probability,
}

#[derive(Clone, Copy, ValueEnum)]
enum Alphabet {
    #[value(name = "ACGT")]
    Acgt,
    #[value(name = "ACGU")]
    Acgu,
    #[value(name = "PROTEIN")]
    Protein,
}

impl Alphabet {
    fn as_str(self) -> &'static str {
        match self {
            Alphabet::Acgt => "ACGT",
            Alphabet::Acgu => "ACGU",
            Alphabet::Protein => "PROTEIN",
        }
    }
}

impl MatrixTypeArg {
    fn to_matrix_type(self) -> MatrixType {
        match self {
            MatrixTypeArg::Auto => MatrixType::Auto,
            MatrixTypeArg::Logodds => MatrixType::Logodds,
            MatrixTypeArg::Probability => MatrixType::Probability,
        }
    }
}

fn parse_positive(value: &str) -> Result<f64, String> {
    let v = value
        .parse::<f64>()
        .map_err(|_| format!("invalid -a value: {}", value))?;
    if v <= 0.0 {
        return Err(format!("-a must be > 0, got {}", v));
    }
    Ok(v)
}

fn parse_background(value: &str) -> Result<Vec<f64>, String> {
    let mut values = Vec::new();
    for part in value.split(',') {
        let v = part
            .parse::<f64>()
            .map_err(|_| format!("invalid -b value: {}", value))?;
        if v <= 0.0 || v > 1.0 {
            return Err(format!("-b values must be in (0, 1], got {}", v));
        }
        values.push(v);
    }
    if values.is_empty() {
        return Err("-b must contain at least one value".to_string());
    }
    if values.len() > 1 {
        let total: f64 = values.iter().sum();
        if (total - 1.0).abs() > 1e-3 {
            return Err(format!("-b vector must sum to 1.0, got {:.6}", total));
        }
    }
    Ok(values)
}

fn parse_nonneg(value: &str) -> Result<f64, String> {
    let v = value
        .parse::<f64>()
        .map_err(|_| format!("invalid value: {}", value))?;
    if v < 0.0 {
        return Err(format!("value must be >= 0, got {}", v));
    }
    Ok(v)
}

fn parse_positive_usize(value: &str) -> Result<usize, String> {
    let v = value
        .parse::<usize>()
        .map_err(|_| format!("invalid --nsites value: {}", value))?;
    if v == 0 {
        return Err("--nsites must be > 0".to_string());
    }
    Ok(v)
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let background = parse_background(&args.background)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let reader: Box<dyn BufRead> = if args.input == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else if args.input.ends_with(".gz") {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(file))
    };

    let matrix_type = args.matrix_type.to_matrix_type();
    let raw_motifs = match args.input_format {
        InputFormat::Json => read_json_with_strict(
            reader,
            args.pseudocount,
            matrix_type,
            &background,
            args.strict,
        )?,
        InputFormat::Homer => read_homer_with_strict(
            reader,
            args.pseudocount,
            matrix_type,
            args.alphabet.as_str(),
            &background,
            args.strict,
        )?,
    };

    let mut processed_motifs = Vec::new();
    for mut m in raw_motifs {
        if !args.extract.is_empty() && m.id != args.extract && m.description != args.extract {
            continue;
        }

        if args.rc {
            if let Err(e) = m.reverse_complement() {
                eprintln!("Warning: skipping motif '{}': {}", m.id, e);
                continue;
            }
        }
        if args.trim_edges > 0.0 {
            m.trim_edges(args.trim_edges);
        }
        if m.matrix.is_empty() {
            continue;
        }
        if args.min_ic > 0.0 && m.total_ic() < args.min_ic {
            continue;
        }
        processed_motifs.push(m);
    }

    let mut stdout = BufWriter::new(io::stdout().lock());
    let output_background = if background.len() > 1 {
        Some(background.as_slice())
    } else {
        None
    };
    write_meme_with_background(
        &mut stdout,
        &processed_motifs,
        args.nsites,
        args.evalue,
        args.renormalize,
        output_background,
        args.strict,
    )?;

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
