//! meme2homer - Convert MEME motif format to HOMER motif format
//!
//! Build (Cargo, recommended - includes .gz support):
//!   cd rust_scripts && cargo build --release
//!   # Binaries: rust_scripts/target/release/meme2homer  homer2meme
//!
//! Install to PATH:
//!   cd rust_scripts && cargo install --path .
//!
//! Usage:
//!   meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
//!   meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer
//!   meme2homer -i motifs.meme -f json > motifs.json
//!   zcat motifs.meme.gz | meme2homer -i -

use clap::{ArgAction, Parser, ValueEnum};
use flate2::read::MultiGzDecoder;
use motif_bridge::io::{read_meme, write_homer, write_json};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};

#[derive(Parser)]
#[command(
    name = "meme2homer",
    about = "Convert MEME format to HOMER motif format.",
    version,
    after_help = "Examples:\n  meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer\n  meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer\n  meme2homer -i motifs.meme -b 0.25 -t 6\n  meme2homer -i motifs.meme -f json > motifs.json\n"
)]
struct Args {
    /// Input MEME file (.meme or .meme.gz), or '-' for stdin
    #[arg(short = 'i', value_name = "<file>")]
    input: String,
    /// Database name (default: NA)
    #[arg(short = 'j', default_value = "NA")]
    db: String,
    /// Override motif name
    #[arg(short = 'k', default_value = "")]
    motif_name: String,
    /// Extract only specified motif by id or name
    #[arg(short = 'e', default_value = "")]
    extract: String,
    /// Background probability scalar or comma-separated vector
    #[arg(short = 'b', default_value = "0.25")]
    bg: String,
    /// Threshold offset in log2 bits
    #[arg(short = 't', default_value = "4.0")]
    t_offset: f64,
    /// Output format: homer (default) or json
    #[arg(short = 'f', long = "format", value_enum, default_value_t = OutputFormat::Homer)]
    output_format: OutputFormat,
    /// Alphabet override for MEME input: ACGT, ACGU, or PROTEIN
    #[arg(long = "alphabet", value_enum)]
    alphabet: Option<Alphabet>,
    /// Output the reverse complement of the motif (DNA/RNA only)
    #[arg(long = "rc", action = ArgAction::SetTrue)]
    rc: bool,
    /// Trim edges with information content below threshold
    #[arg(long = "trim-edges", value_parser = parse_nonneg, default_value = "0.0")]
    trim_edges: f64,
    /// Filter out motifs with total information content below threshold
    #[arg(long = "min-ic", value_parser = parse_nonneg, default_value = "0.0")]
    min_ic: f64,
    /// Renormalize each row before writing HOMER output
    #[arg(long = "renormalize", action = ArgAction::SetTrue)]
    renormalize: bool,
    /// Keep an existing motif threshold when present
    #[arg(long = "keep-threshold", action = ArgAction::SetTrue)]
    keep_threshold: bool,
}

#[derive(Clone, Copy, ValueEnum)]
enum OutputFormat {
    #[value(name = "homer")]
    Homer,
    #[value(name = "json")]
    Json,
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

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let bg =
        parse_background(&args.bg).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let reader: Box<dyn BufRead> = if args.input == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else if args.input.ends_with(".gz") {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(file))
    };

    let raw_motifs = read_meme(reader, args.alphabet.map(|a| a.as_str()))?;
    let mut processed_motifs = Vec::new();

    for mut m in raw_motifs {
        let original_name = m.description.clone();
        if !args.extract.is_empty() && m.id != args.extract && original_name != args.extract {
            continue;
        }

        m.description = if !args.motif_name.is_empty() {
            format!("{}/{}", args.motif_name, args.db)
        } else {
            format!("{}/{}", original_name, args.db)
        };

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
    match args.output_format {
        OutputFormat::Json => write_json(&mut stdout, &processed_motifs)?,
        OutputFormat::Homer => write_homer(
            &mut stdout,
            &processed_motifs,
            &bg,
            args.t_offset,
            args.keep_threshold,
            args.renormalize,
        )?,
    }

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
