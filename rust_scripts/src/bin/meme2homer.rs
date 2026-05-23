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

use flate2::read::MultiGzDecoder;
use motif_bridge::io::{read_meme, write_homer, write_json};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};

// ---------------------------------------------------------------------------
// Argument parsing
// ---------------------------------------------------------------------------

struct Args {
    input: String,
    db: String,
    motif_name: String,
    extract: String,
    bg: f64,
    t_offset: f64,
    output_format: OutputFormat,
    alphabet: String,
    rc: bool,
    trim_edges: f64,
    min_ic: f64,
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Homer,
    Json,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut input = String::new();
    let mut db = String::from("NA");
    let mut motif_name = String::new();
    let mut extract = String::new();
    let mut bg = 0.25_f64;
    let mut t_offset = 4.0_f64;
    let mut output_format = OutputFormat::Homer;
    let mut alphabet = String::from("ACGT");
    let mut rc = false;
    let mut trim_edges = 0.0_f64;
    let mut min_ic = 0.0_f64;
    let mut i = 1usize;

    while i < argv.len() {
        match argv[i].as_str() {
            "-i" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-i requires an argument.".into());
                }
                input = argv[i].clone();
            }
            "-j" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-j requires an argument.".into());
                }
                db = argv[i].clone();
            }
            "-k" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-k requires an argument.".into());
                }
                motif_name = argv[i].clone();
            }
            "-e" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-e requires an argument.".into());
                }
                extract = argv[i].clone();
            }
            "-b" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-b requires an argument.".into());
                }
                bg = argv[i]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid -b value: {}", argv[i]))?;
                if bg <= 0.0 || bg > 1.0 {
                    return Err(format!("-b must be in (0, 1], got {}", bg));
                }
            }
            "-t" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-t requires an argument.".into());
                }
                t_offset = argv[i]
                    .parse::<f64>()
                    .map_err(|_| format!("invalid -t value: {}", argv[i]))?;
            }
            "-f" | "--format" => {
                i += 1;
                if i >= argv.len() {
                    return Err(format!("{} requires an argument.", argv[i - 1]));
                }
                output_format = match argv[i].as_str() {
                    "homer" => OutputFormat::Homer,
                    "json" => OutputFormat::Json,
                    other => return Err(format!("unknown format: {}", other)),
                };
            }
            "--alphabet" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--alphabet requires an argument.".into());
                }
                let val = argv[i].as_str();
                if val != "ACGT" && val != "ACGU" && val != "PROTEIN" {
                    return Err(format!("unknown alphabet: {}", val));
                }
                alphabet = val.to_string();
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
        db,
        motif_name,
        extract,
        bg,
        t_offset,
        output_format,
        alphabet,
        rc,
        trim_edges,
        min_ic,
    })
}

fn usage() {
    eprintln!(
        r#"Usage: meme2homer -i <input_file> [OPTIONS]

Convert MEME format to HOMER motif format.

Options:
    -i <file>    Input MEME file (.meme or .meme.gz), or '-' for stdin
    -j <string>  Database name (default: NA)
    -k <string>  Override motif name
    -e <string>  Extract only specified motif by id or name
    -b <float>   Background probability in (0, 1] (default: 0.25)
    -t <float>   Threshold offset in log2 bits (default: 4.0)
    -f, --format <fmt>  Output format: homer (default) or json
    --alphabet <str> Alphabet: ACGT (DNA, default), ACGU (RNA), or PROTEIN
    --rc         Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float> Filter out motifs with total information content below threshold
    -h           Show this help

Examples:
    meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme -b 0.25 -t 6
    meme2homer -i motifs.meme -f json > motifs.json
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

    let raw_motifs = read_meme(reader, &args.alphabet)?;
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
        OutputFormat::Homer => write_homer(&mut stdout, &processed_motifs, args.bg, args.t_offset)?,
    }

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
