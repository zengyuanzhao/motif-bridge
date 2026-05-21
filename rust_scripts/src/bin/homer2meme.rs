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
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Clone, Copy)]
enum InputFormat {
    Homer,
    Json,
}

#[derive(Clone, Copy)]
enum MatrixType {
    Auto,
    Logodds,
    Probability,
}

struct Args {
    input: String,
    extract: String,
    pseudocount: f64,
    input_format: InputFormat,
    matrix_type: MatrixType,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut input = String::new();
    let mut extract = String::new();
    let mut pseudocount = 0.01_f64;
    let mut input_format = InputFormat::Homer;
    let mut matrix_type = MatrixType::Auto;
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
            "-f" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-f requires a value".into());
                }
                input_format = match argv[i].as_str() {
                    "homer" => InputFormat::Homer,
                    "json" => InputFormat::Json,
                    other => return Err(format!("unknown format: {}", other)),
                };
            }
            "--format" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--format requires a value".into());
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

fn is_logodds(row: &[f64], matrix_type: MatrixType) -> bool {
    match matrix_type {
        MatrixType::Logodds => true,
        MatrixType::Probability => false,
        MatrixType::Auto => {
            let s: f64 = row.iter().sum();
            !(0.98..=1.02).contains(&s)
        }
    }
}

fn logodds_to_prob(row: &[f64], pseudocount: f64) -> Vec<f64> {
    let background = 0.25_f64;
    let raw: Vec<f64> = row.iter().map(|&v| 2.0_f64.powf(v) * background).collect();
    let total: f64 = raw.iter().sum::<f64>() + pseudocount * raw.len() as f64;
    raw.iter().map(|&v| (v + pseudocount) / total).collect()
}

fn print_meme_header() {
    println!("MEME version 4");
    println!();
    println!("ALPHABET= ACGT");
    println!();
    println!("strands: + -");
    println!();
    println!("Background letter frequencies");
    println!("A 0.25 C 0.25 G 0.25 T 0.25");
    println!();
}

fn print_meme_motif(id: &str, desc: &str, matrix: &[Vec<f64>]) {
    let width = matrix.len();
    println!("MOTIF {} {}", id, desc);
    println!();
    println!(
        "letter-probability matrix: alength= 4 w= {} nsites= 20 E= 0",
        width
    );
    for row in matrix {
        let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
        println!("  {}", line.join("  "));
    }
    println!();
}

fn parse_and_convert_homer<R: BufRead>(reader: R, args: &Args) -> io::Result<()> {
    let mut header_printed = false;
    let mut in_motif = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if in_motif && !matrix.is_empty() {
                if !header_printed {
                    print_meme_header();
                    header_printed = true;
                }
                print_meme_motif(&motif_id, &description, &matrix);
            }
            matrix.clear();

            let parts: Vec<&str> = rest.splitn(6, '\t').collect();
            let mid = parts.first().copied().unwrap_or("motif").to_string();
            let raw_desc = parts.get(1).copied().unwrap_or("").to_string();
            let final_desc = if raw_desc.is_empty() {
                mid.clone()
            } else {
                raw_desc
            };

            if !args.extract.is_empty() && mid != args.extract && final_desc != args.extract {
                in_motif = false;
                continue;
            }
            motif_id = mid;
            description = final_desc;
            in_motif = true;
            continue;
        }

        if !in_motif {
            continue;
        }

        let row_result: Result<Vec<f64>, _> = trimmed
            .split_whitespace()
            .map(|s| s.parse::<f64>())
            .collect();
        if let Ok(mut row) = row_result {
            if !row.is_empty() {
                if row.len() != 4 {
                    eprintln!(
                        "Warning: skipping malformed matrix row with {} cols (expected 4): {}",
                        row.len(),
                        trimmed
                    );
                    continue;
                }
                if is_logodds(&row, args.matrix_type) {
                    row = logodds_to_prob(&row, args.pseudocount);
                }
                matrix.push(row);
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        if !header_printed {
            print_meme_header();
        }
        print_meme_motif(&motif_id, &description, &matrix);
    }

    Ok(())
}

fn parse_and_convert_json<R: BufRead>(reader: R, args: &Args) -> io::Result<()> {
    let mut content = String::new();
    for line_result in reader.lines() {
        content.push_str(&line_result?);
    }

    let data: serde_json::Value = match serde_json::from_str(&content) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Error: invalid JSON: {}", e);
            std::process::exit(1);
        }
    };

    let mut header_printed = false;

    if let Some(motifs) = data.get("motifs").and_then(|m| m.as_array()) {
        for motif in motifs {
            let mid = motif.get("id").and_then(|v| v.as_str()).unwrap_or("motif");
            let desc = motif
                .get("description")
                .and_then(|v| v.as_str())
                .unwrap_or(mid);
            let matrix = motif.get("matrix").and_then(|v| v.as_array());

            if matrix.is_none() {
                continue;
            }
            let matrix = matrix.unwrap();

            if !args.extract.is_empty() && mid != args.extract && desc != args.extract {
                continue;
            }

            let mut processed: Vec<Vec<f64>> = Vec::new();
            for row in matrix {
                if let Some(arr) = row.as_array() {
                    if arr.len() != 4 {
                        eprintln!("Warning: skipping malformed matrix row (expected 4 cols)");
                        continue;
                    }
                    let vals: Vec<f64> = arr.iter().filter_map(|v| v.as_f64()).collect();
                    if vals.len() == 4 {
                        if is_logodds(&vals, MatrixType::Auto) {
                            processed.push(logodds_to_prob(&vals, args.pseudocount));
                        } else {
                            processed.push(vals);
                        }
                    }
                }
            }

            if !processed.is_empty() {
                if !header_printed {
                    print_meme_header();
                    header_printed = true;
                }
                print_meme_motif(mid, desc, &processed);
            }
        }
    }

    Ok(())
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

    match args.input_format {
        InputFormat::Json => parse_and_convert_json(reader, &args)?,
        InputFormat::Homer => parse_and_convert_homer(reader, &args)?,
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
