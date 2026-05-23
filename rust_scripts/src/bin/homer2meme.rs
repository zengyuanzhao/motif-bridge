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

use flate2::read::MultiGzDecoder;
use motif_bridge::io::{read_homer, read_json, write_meme, MatrixType};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};

#[derive(Clone, Copy)]
enum InputFormat {
    Homer,
    Json,
}

struct Args {
    input: String,
    extract: String,
    pseudocount: f64,
    input_format: InputFormat,
    matrix_type: MatrixType,
    rc: bool,
    trim_edges: f64,
    min_ic: f64,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut input = String::new();
    let mut extract = String::new();
    let mut pseudocount = 0.01_f64;
    let mut input_format = InputFormat::Homer;
    let mut matrix_type = MatrixType::Auto;
    let mut rc = false;
    let mut trim_edges = 0.0_f64;
    let mut min_ic = 0.0_f64;
    let mut i = 1usize;

    while i < argv.len() {
        match argv[i].as_str() {
            "-i" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-i requires a value".into());
                }
                input = argv[i].clone();
            }
            "-e" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-e requires a value".into());
                }
                extract = argv[i].clone();
            }
            "-a" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-a requires a value".into());
                }
                pseudocount = argv[i]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid -a value: {}", argv[i]))?;
                if pseudocount <= 0.0 {
                    return Err(format!("-a must be > 0, got {}", pseudocount));
                }
            }
            "-f" | "--format" => {
                i += 1;
                if i >= argv.len() {
                    return Err(format!("{} requires a value", argv[i - 1]));
                }
                input_format = match argv[i].as_str() {
                    "homer" => InputFormat::Homer,
                    "json" => InputFormat::Json,
                    other => return Err(format!("unknown format: {}", other)),
                };
            }
            "--input-format" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--input-format requires a value".into());
                }
                matrix_type = match argv[i].as_str() {
                    "auto" => MatrixType::Auto,
                    "logodds" => MatrixType::Logodds,
                    "probability" => MatrixType::Probability,
                    other => return Err(format!("unknown input-format: {}", other)),
                };
            }
            "--rc" => {
                rc = true;
            }
            "--trim-edges" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--trim-edges requires an argument.".into());
                }
                trim_edges = argv[i]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid --trim-edges value: {}", argv[i]))?;
                if trim_edges < 0.0 {
                    return Err(format!("--trim-edges must be >= 0, got {}", trim_edges));
                }
            }
            "--min-ic" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--min-ic requires an argument.".into());
                }
                min_ic = argv[i]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid --min-ic value: {}", argv[i]))?;
                if min_ic < 0.0 {
                    return Err(format!("--min-ic must be >= 0, got {}", min_ic));
                }
            }
            "-h" | "--help" => {
                usage();
                std::process::exit(0);
            }
            other => {
                eprintln!("Unknown option: {}", other);
                usage();
                std::process::exit(1);
            }
        }
        i += 1;
    }
    if input.is_empty() {
        return Err("-i <input_file> is required.".into());
    }
    Ok(Args {
        input,
        extract,
        pseudocount,
        input_format,
        matrix_type,
        rc,
        trim_edges,
        min_ic,
    })
}

fn usage() {
    eprintln!(
        r#"Usage: homer2meme -i <input_file> [OPTIONS]

Convert HOMER motif format to MEME format.

Options:
    -i <file>    Input HOMER motif file (.homer, .motif, .json, or .gz), or '-' for stdin
    -e <string>  Extract only specified motif by id or description
    -a <float>   Pseudocount for log-odds -> probability conversion (default: 0.01)
    -f, --format <fmt>  Input format: homer (default) or json
    --input-format <fmt>  Matrix type: auto (default), logodds, or probability
    --rc         Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float> Filter out motifs with total information content below threshold
    -h           Show this help

Examples:
    homer2meme -i motifs.homer > motifs.meme
    homer2meme -i motifs.homer.gz > motifs.meme
    homer2meme -i motifs.homer -e "CTCF/Jaspar"
    homer2meme -i motifs.json -f json > motifs.meme
    homer2meme -i motifs.homer --input-format logodds
    cat motifs.homer | homer2meme -i -
"#
    );
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let argv: Vec<String> = env::args().collect();
    let args = match parse_args(&argv) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            usage();
            std::process::exit(1);
        }
    };

    let reader: Box<dyn BufRead> = if args.input == "-" {
        Box::new(BufReader::new(io::stdin().lock()))
    } else if args.input.ends_with(".gz") {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        let file = File::open(&args.input)?;
        Box::new(BufReader::new(file))
    };

    let raw_motifs = match args.input_format {
        InputFormat::Json => read_json(reader, args.pseudocount)?,
        InputFormat::Homer => read_homer(reader, args.pseudocount, args.matrix_type)?,
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
    write_meme(&mut stdout, &processed_motifs)?;

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
