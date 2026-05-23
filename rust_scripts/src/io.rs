use crate::{Motif, MotifError};
use serde_json::Value;
use std::io::{BufRead, Write};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MatrixType {
    Auto,
    Logodds,
    Probability,
}

pub fn is_logodds(row: &[f64], matrix_type: MatrixType) -> bool {
    match matrix_type {
        MatrixType::Logodds => true,
        MatrixType::Probability => false,
        MatrixType::Auto => {
            let s: f64 = row.iter().sum();
            !(0.98..=1.02).contains(&s)
        }
    }
}

pub fn logodds_to_prob(row: &[f64], pseudocount: f64) -> Vec<f64> {
    let background = 0.25;
    let raw: Vec<f64> = row.iter().map(|&v| 2.0_f64.powf(v) * background).collect();
    let total: f64 = raw.iter().sum::<f64>() + pseudocount * raw.len() as f64;
    raw.into_iter().map(|v| (v + pseudocount) / total).collect()
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

pub fn read_meme<R: BufRead>(reader: R, alphabet_arg: &str) -> Result<Vec<Motif>, MotifError> {
    let mut in_motif = false;
    let mut in_matrix = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut motifs = Vec::new();

    let mut expected_cols = match alphabet_arg {
        "ACGT" | "ACGU" => 4,
        "PROTEIN" => 20,
        other => other.len(),
    };

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("MOTIF") {
            if in_motif && !matrix.is_empty() {
                motifs.push(Motif::new(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    alphabet_arg.to_string(),
                ));
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

            motif_id = id;
            description = original_name;
            matrix.clear();
            in_motif = true;
            continue;
        }

        if !in_motif {
            continue;
        }

        if trimmed.starts_with("URL") {
            continue;
        }

        if trimmed.starts_with("//") {
            if !matrix.is_empty() {
                motifs.push(Motif::new(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    alphabet_arg.to_string(),
                ));
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
                }
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        motifs.push(Motif::new(
            motif_id,
            description,
            matrix,
            alphabet_arg.to_string(),
        ));
    }

    Ok(motifs)
}

pub fn read_homer<R: BufRead>(
    reader: R,
    pseudocount: f64,
    matrix_type: MatrixType,
) -> Result<Vec<Motif>, MotifError> {
    let mut in_motif = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut motifs = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if in_motif && !matrix.is_empty() {
                motifs.push(Motif::new(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    "ACGT".to_string(),
                ));
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
                    continue; // Skip malformed
                }
                if is_logodds(&row, matrix_type) {
                    row = logodds_to_prob(&row, pseudocount);
                }
                matrix.push(row);
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        motifs.push(Motif::new(
            motif_id,
            description,
            matrix,
            "ACGT".to_string(),
        ));
    }

    Ok(motifs)
}

pub fn read_json<R: BufRead>(mut reader: R, pseudocount: f64) -> Result<Vec<Motif>, MotifError> {
    let mut content = String::new();
    reader.read_to_string(&mut content)?;

    let data: Value = serde_json::from_str(&content)?;
    let mut motifs = Vec::new();

    if let Some(motifs_array) = data.get("motifs").and_then(|m| m.as_array()) {
        for motif in motifs_array {
            let mid = motif.get("id").and_then(|v| v.as_str()).unwrap_or("motif");
            let desc = motif
                .get("description")
                .and_then(|v| v.as_str())
                .unwrap_or(mid);
            let matrix = motif.get("matrix").and_then(|v| v.as_array());
            let alphabet = motif
                .get("alphabet")
                .and_then(|v| v.as_str())
                .unwrap_or("ACGT");

            if matrix.is_none() {
                continue;
            }
            let matrix = matrix.unwrap();

            let mut processed: Vec<Vec<f64>> = Vec::new();
            for row in matrix {
                if let Some(arr) = row.as_array() {
                    if arr.len() != 4 {
                        continue;
                    }
                    let vals: Vec<f64> = arr.iter().filter_map(|v| v.as_f64()).collect();
                    if vals.len() == 4 {
                        if is_logodds(&vals, MatrixType::Auto) {
                            processed.push(logodds_to_prob(&vals, pseudocount));
                        } else {
                            processed.push(vals);
                        }
                    }
                }
            }

            if !processed.is_empty() {
                motifs.push(Motif::new(
                    mid.to_string(),
                    desc.to_string(),
                    processed,
                    alphabet.to_string(),
                ));
            }
        }
    }

    Ok(motifs)
}

pub fn write_homer<W: Write>(
    writer: &mut W,
    motifs: &[Motif],
    bg: f64,
    t_offset: f64,
) -> Result<(), MotifError> {
    for motif in motifs {
        let score = motif.calculate_score(bg, t_offset);
        writeln!(
            writer,
            ">{}\t{}\t{:.6}\t0\t0\t0",
            motif.id, motif.description, score
        )?;
        for row in &motif.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            writeln!(writer, "{}", line.join("\t"))?;
        }
    }
    Ok(())
}

pub fn write_meme<W: Write>(writer: &mut W, motifs: &[Motif]) -> Result<(), MotifError> {
    let mut header_printed = false;
    for motif in motifs {
        if !header_printed {
            writeln!(writer, "MEME version 4\n")?;
            writeln!(writer, "ALPHABET= {}\n", motif.alphabet)?;
            writeln!(writer, "strands: + -\n")?;
            writeln!(writer, "Background letter frequencies")?;
            writeln!(writer, "A 0.25 C 0.25 G 0.25 T 0.25\n")?;
            header_printed = true;
        }

        let expected_cols = if motif.alphabet == "PROTEIN" { 20 } else { 4 };
        let width = motif.matrix.len();
        writeln!(writer, "MOTIF {} {}", motif.id, motif.description)?;
        writeln!(writer)?;
        writeln!(
            writer,
            "letter-probability matrix: alength= {} w= {} nsites= 20 E= 0",
            expected_cols, width
        )?;
        for row in &motif.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            writeln!(writer, "  {}", line.join("  "))?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

pub fn write_json<W: Write>(writer: &mut W, motifs: &[Motif]) -> Result<(), MotifError> {
    writeln!(writer, "{{")?;
    writeln!(writer, "  \"version\": \"1.0\",")?;
    writeln!(writer, "  \"source\": \"meme\",")?;
    writeln!(writer, "  \"motifs\": [")?;
    for (mi, motif) in motifs.iter().enumerate() {
        writeln!(writer, "    {{")?;
        writeln!(writer, "      \"id\": \"{}\",", escape_json(&motif.id))?;
        writeln!(
            writer,
            "      \"description\": \"{}\",",
            escape_json(&motif.description)
        )?;
        if !motif.alphabet.is_empty() {
            writeln!(
                writer,
                "      \"alphabet\": \"{}\",",
                escape_json(&motif.alphabet)
            )?;
        }
        writeln!(writer, "      \"matrix\": [")?;
        for (ri, row) in motif.matrix.iter().enumerate() {
            let vals: Vec<String> = row.iter().map(|v| format!("{}", v)).collect();
            if ri + 1 < motif.matrix.len() {
                writeln!(writer, "        [{}],", vals.join(", "))?;
            } else {
                writeln!(writer, "        [{}]", vals.join(", "))?;
            }
        }
        writeln!(writer, "      ]")?;
        if mi + 1 < motifs.len() {
            writeln!(writer, "    }},")?;
        } else {
            writeln!(writer, "    }}")?;
        }
    }
    writeln!(writer, "  ]")?;
    writeln!(writer, "}}")?;
    Ok(())
}
