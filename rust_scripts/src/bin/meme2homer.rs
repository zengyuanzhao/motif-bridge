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
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

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
            "-f" => {
                i += 1;
                if i >= argv.len() {
                    return Err("-f requires an argument.".into());
                }
                output_format = match argv[i].as_str() {
                    "homer" => OutputFormat::Homer,
                    "json" => OutputFormat::Json,
                    other => return Err(format!("unknown format: {}", other)),
                };
            }
            "--format" => {
                i += 1;
                if i >= argv.len() {
                    return Err("--format requires an argument.".into());
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
    -h           Show this help

Examples:
    meme2homer -i motifs.meme -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme.gz -j JASPAR2026 > motifs.homer
    meme2homer -i motifs.meme -b 0.25 -t 6
    meme2homer -i motifs.meme -f json > motifs.json
"#
    );
}

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct Motif {
    id: String,
    description: String,
    matrix: Vec<Vec<f64>>,
    alphabet: String,
}

impl Motif {
    fn calculate_score(&self, bg: f64, t_offset: f64) -> f64 {
        let raw: f64 = self
            .matrix
            .iter()
            .map(|row| {
                let max_p = row.iter().cloned().fold(0.0_f64, f64::max);
                if max_p > 0.0 {
                    (max_p / bg).log2()
                } else {
                    0.0
                }
            })
            .sum();
        (raw - t_offset).max(0.0)
    }

    fn print_homer(&self, bg: f64, t_offset: f64) {
        let score = self.calculate_score(bg, t_offset);
        println!(">{}\t{}\t{:.6}\t0\t0\t0", self.id, self.description, score);
        for row in &self.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            println!("{}", line.join("\t"));
        }
    }
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

fn print_json(motifs: &[Motif]) {
    println!("{{");
    println!("  \"version\": \"1.0\",");
    println!("  \"source\": \"meme\",");
    println!("  \"motifs\": [");
    for (mi, motif) in motifs.iter().enumerate() {
        println!("    {{");
        println!("      \"id\": \"{}\",", escape_json(&motif.id));
        println!(
            "      \"description\": \"{}\",",
            escape_json(&motif.description)
        );
        if !motif.alphabet.is_empty() {
            println!("      \"alphabet\": \"{}\",", escape_json(&motif.alphabet));
        }
        println!("      \"matrix\": [");
        for (ri, row) in motif.matrix.iter().enumerate() {
            let vals: Vec<String> = row.iter().map(|v| format!("{}", v)).collect();
            if ri + 1 < motif.matrix.len() {
                println!("        [{}],", vals.join(", "));
            } else {
                println!("        [{}]", vals.join(", "));
            }
        }
        println!("      ]");
        if mi + 1 < motifs.len() {
            println!("    }},");
        } else {
            println!("    }}");
        }
    }
    println!("  ]");
    println!("}}");
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

// ---------------------------------------------------------------------------
// Core parsing
// ---------------------------------------------------------------------------

fn parse_and_convert<R: BufRead>(reader: R, args: &Args) -> io::Result<()> {
    let mut in_motif = false;
    let mut in_matrix = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut motifs: Vec<Motif> = Vec::new();

    let mut expected_cols = match args.alphabet.as_str() {
        "ACGT" | "ACGU" => 4,
        "PROTEIN" => 20,
        other => other.len(),
    };

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("MOTIF") {
            if in_motif && !matrix.is_empty() {
                let m = Motif {
                    id: motif_id.clone(),
                    description: description.clone(),
                    matrix: std::mem::take(&mut matrix),
                    alphabet: args.alphabet.clone(),
                };
                if matches!(args.output_format, OutputFormat::Json) {
                    motifs.push(m);
                } else {
                    m.print_homer(args.bg, args.t_offset);
                }
            }
            in_matrix = false;

            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.is_empty() {
                in_motif = false;
                matrix.clear();
                continue;
            }
            let id = parts[0].to_string();
            let original_name = if parts.len() > 1 {
                parts[1..].join(" ")
            } else {
                id.clone()
            };

            if !args.extract.is_empty() && id != args.extract && original_name != args.extract {
                in_motif = false;
                continue;
            }
            motif_id = id;
            description = if !args.motif_name.is_empty() {
                format!("{}/{}", args.motif_name, args.db)
            } else {
                format!("{}/{}", original_name, args.db)
            };
            in_motif = true;
            continue;
        }

        if !in_motif {
            continue;
        }

        if trimmed.starts_with("URL ") || trimmed == "URL" {
            continue;
        }

        if trimmed.starts_with("//") {
            if !matrix.is_empty() {
                let m = Motif {
                    id: motif_id.clone(),
                    description: description.clone(),
                    matrix: std::mem::take(&mut matrix),
                    alphabet: args.alphabet.clone(),
                };
                if matches!(args.output_format, OutputFormat::Json) {
                    motifs.push(m);
                } else {
                    m.print_homer(args.bg, args.t_offset);
                }
            }
            in_motif = false;
            in_matrix = false;
            continue;
        }

        if trimmed.starts_with("letter-probability matrix:") {
            in_matrix = true;
            if let Some(alength_idx) = trimmed.find("alength=") {
                let rest = &trimmed[alength_idx + 8..];
                if let Some(space_idx) = rest.find(char::is_whitespace) {
                    if let Ok(len) = rest[..space_idx].parse::<usize>() {
                        expected_cols = len;
                    }
                } else if let Ok(len) = rest.parse::<usize>() {
                    expected_cols = len;
                }
            }
            continue;
        }

        if in_matrix {
            let first = trimmed.chars().next();
            if !matches!(first, Some('0'..='9') | Some('.')) {
                continue;
            }
            let row: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();
            if let Ok(values) = row {
                if values.len() == expected_cols {
                    matrix.push(values);
                } else if !values.is_empty() {
                    eprintln!(
                        "Warning: skipping malformed matrix row with {} cols (expected {})",
                        values.len(),
                        expected_cols
                    );
                }
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        let m = Motif {
            id: motif_id,
            description,
            matrix,
            alphabet: args.alphabet.clone(),
        };
        if matches!(args.output_format, OutputFormat::Json) {
            motifs.push(m);
        } else {
            m.print_homer(args.bg, args.t_offset);
        }
    }

    if matches!(args.output_format, OutputFormat::Json) {
        print_json(&motifs);
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

    parse_and_convert(reader, &args)?;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
