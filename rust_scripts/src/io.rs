use crate::{background_values, renormalized, Motif, MotifError};
use serde_json::Value;
use std::io::{BufRead, Write};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MatrixType {
    Auto,
    Logodds,
    Probability,
}

fn alphabet_letters(alphabet: &str) -> &str {
    match alphabet {
        "ACGT" => "ACGT",
        "ACGU" => "ACGU",
        "PROTEIN" => "ACDEFGHIKLMNPQRSTVWY",
        other => other,
    }
}

fn background_line(alphabet: &str) -> String {
    let letters = alphabet_letters(alphabet);
    let count = letters.chars().count();
    if count == 0 {
        return String::new();
    }
    let freq = 1.0_f64 / count as f64;
    let mut parts = Vec::with_capacity(count);
    for ch in letters.chars() {
        let formatted = if count == 4 || count == 20 {
            format!("{:.2}", freq)
        } else {
            format!("{:.6}", freq)
        };
        parts.push(format!("{ch} {formatted}"));
    }
    parts.join(" ")
}

fn parse_header_value<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let idx = line.find(key)?;
    line[idx + key.len()..].split_whitespace().next()
}

pub fn is_logodds(row: &[f64], matrix_type: MatrixType) -> bool {
    match matrix_type {
        MatrixType::Logodds => true,
        MatrixType::Probability => false,
        MatrixType::Auto => {
            let s: f64 = row.iter().sum();
            let is_probability = (0.98..=1.02).contains(&s);
            if !is_probability && row.iter().all(|v| *v >= 0.0) && (0.90..=1.10).contains(&s) {
                eprintln!(
                    "Warning: row sum {:.4} is close to the auto-detection boundary; use --input-format logodds or --input-format probability for reproducibility",
                    s
                );
            }
            !is_probability
        }
    }
}

pub fn logodds_to_prob(
    row: &[f64],
    pseudocount: f64,
    background: &[f64],
) -> Result<Vec<f64>, MotifError> {
    let bg_values = background_values(background, row.len())?;
    let raw: Vec<f64> = row
        .iter()
        .enumerate()
        .map(|(i, &v)| 2.0_f64.powf(v) * bg_values[i])
        .collect();
    let total: f64 = raw.iter().sum::<f64>() + pseudocount * raw.len() as f64;
    Ok(raw.into_iter().map(|v| (v + pseudocount) / total).collect())
}

fn json_string(value: &str) -> String {
    serde_json::to_string(value).unwrap_or_else(|_| format!("{:?}", value))
}

pub fn read_meme<R: BufRead>(
    reader: R,
    alphabet_override: Option<&str>,
) -> Result<Vec<Motif>, MotifError> {
    let mut in_motif = false;
    let mut in_matrix = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut nsites: Option<usize> = None;
    let mut evalue: Option<f64> = None;
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut motifs = Vec::new();

    let override_set = alphabet_override.is_some();
    let mut alphabet = alphabet_override.unwrap_or("ACGT").to_string();
    let mut expected_cols = alphabet_letters(&alphabet).chars().count();

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("ALPHABET=") {
            if !override_set {
                let detected = rest.split_whitespace().next().unwrap_or("");
                if !detected.is_empty() {
                    alphabet = detected.to_string();
                    expected_cols = alphabet_letters(&alphabet).chars().count();
                }
            }
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix("MOTIF") {
            if !rest.is_empty()
                && !rest
                    .chars()
                    .next()
                    .map(|c| c.is_whitespace())
                    .unwrap_or(false)
            {
                in_matrix = false;
                continue;
            }
            if in_motif && !matrix.is_empty() {
                motifs.push(Motif::with_metadata(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    alphabet.clone(),
                    None,
                    nsites,
                    evalue,
                ));
            }
            in_matrix = false;

            let rest = rest.trim();
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
            nsites = None;
            evalue = None;
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
                motifs.push(Motif::with_metadata(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    alphabet.clone(),
                    None,
                    nsites,
                    evalue,
                ));
            }
            in_motif = false;
            in_matrix = false;
            continue;
        }

        if trimmed.starts_with("letter-probability matrix:") {
            in_matrix = true;
            if let Some(value) = parse_header_value(trimmed, "nsites=") {
                nsites = value.parse::<usize>().ok();
            }
            if let Some(value) = parse_header_value(trimmed, "E=") {
                evalue = value.parse::<f64>().ok();
            }
            if let Some(alength_idx) = trimmed.find("alength=") {
                let rest = trimmed[alength_idx + 8..].trim_start();
                if let Some(value) = rest.split_whitespace().next() {
                    let Ok(len) = value.parse::<usize>() else {
                        continue;
                    };
                    if len != expected_cols {
                        eprintln!(
                            "Warning: alength={} conflicts with alphabet {} (expected {} cols); using alphabet-derived width",
                            len, alphabet, expected_cols
                        );
                    } else {
                        expected_cols = len;
                    }
                }
            }
            continue;
        }

        if in_matrix {
            let first = trimmed.chars().next();
            if !matches!(first, Some('0'..='9') | Some('.') | Some('-')) {
                continue;
            }
            let row: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();
            if let Ok(values) = row {
                if values.len() == expected_cols {
                    if values.iter().any(|v| *v < 0.0) {
                        eprintln!(
                            "Warning: negative value in matrix row (expected probabilities): {}",
                            trimmed
                        );
                    }
                    matrix.push(values);
                } else if !values.is_empty() {
                    eprintln!(
                        "Warning: skipping malformed matrix row (expected {} cols, got {}): {}",
                        expected_cols,
                        values.len(),
                        trimmed
                    );
                }
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        motifs.push(Motif::with_metadata(
            motif_id,
            description,
            matrix,
            alphabet,
            None,
            nsites,
            evalue,
        ));
    }

    Ok(motifs)
}

pub fn read_homer<R: BufRead>(
    reader: R,
    pseudocount: f64,
    matrix_type: MatrixType,
    alphabet: &str,
    background: &[f64],
) -> Result<Vec<Motif>, MotifError> {
    let mut in_motif = false;
    let mut motif_id = String::new();
    let mut description = String::new();
    let mut threshold: Option<f64> = None;
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut motifs = Vec::new();
    let expected_cols = alphabet_letters(alphabet).chars().count();

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if in_motif && !matrix.is_empty() {
                motifs.push(Motif::with_threshold(
                    motif_id.clone(),
                    description.clone(),
                    std::mem::take(&mut matrix),
                    alphabet.to_string(),
                    threshold,
                ));
            }
            matrix.clear();

            let parts: Vec<&str> = rest.splitn(6, '\t').collect();
            let mid = parts.first().copied().unwrap_or("motif").to_string();
            let raw_desc = parts.get(1).copied().unwrap_or("").to_string();
            let source_threshold = parts.get(2).and_then(|v| v.parse::<f64>().ok());
            let final_desc = if raw_desc.is_empty() {
                mid.clone()
            } else {
                raw_desc
            };

            motif_id = mid;
            description = final_desc;
            threshold = source_threshold;
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
                if row.len() != expected_cols {
                    eprintln!(
                        "Warning: skipping malformed row (expected {} cols, got {}): {}",
                        expected_cols,
                        row.len(),
                        trimmed
                    );
                    continue; // Skip malformed
                }
                if is_logodds(&row, matrix_type) {
                    row = logodds_to_prob(&row, pseudocount, background)?;
                }
                matrix.push(row);
            }
        }
    }

    if in_motif && !matrix.is_empty() {
        motifs.push(Motif::with_threshold(
            motif_id,
            description,
            matrix,
            alphabet.to_string(),
            threshold,
        ));
    }

    Ok(motifs)
}

pub fn read_json<R: BufRead>(
    mut reader: R,
    pseudocount: f64,
    matrix_type: MatrixType,
    background: &[f64],
) -> Result<Vec<Motif>, MotifError> {
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
            let threshold = motif.get("threshold").and_then(|v| v.as_f64());
            let nsites = motif
                .get("nsites")
                .and_then(|v| v.as_u64())
                .and_then(|v| usize::try_from(v).ok());
            let evalue = motif.get("evalue").and_then(|v| v.as_f64());

            if matrix.is_none() {
                continue;
            }
            let matrix = matrix.unwrap();

            let mut processed: Vec<Vec<f64>> = Vec::new();
            let expected_cols = alphabet_letters(alphabet).chars().count();
            for row in matrix {
                if let Some(arr) = row.as_array() {
                    if arr.len() != expected_cols {
                        continue;
                    }
                    let vals: Vec<f64> = arr.iter().filter_map(|v| v.as_f64()).collect();
                    if vals.len() == expected_cols {
                        if is_logodds(&vals, matrix_type) {
                            processed.push(logodds_to_prob(&vals, pseudocount, background)?);
                        } else {
                            processed.push(vals);
                        }
                    }
                }
            }

            if !processed.is_empty() {
                motifs.push(Motif::with_metadata(
                    mid.to_string(),
                    desc.to_string(),
                    processed,
                    alphabet.to_string(),
                    threshold,
                    nsites,
                    evalue,
                ));
            }
        }
    }

    Ok(motifs)
}

pub fn write_homer<W: Write>(
    writer: &mut W,
    motifs: &[Motif],
    bg: &[f64],
    t_offset: f64,
    keep_threshold: bool,
    renormalize_rows: bool,
) -> Result<(), MotifError> {
    for motif in motifs {
        let score = if keep_threshold {
            if let Some(threshold) = motif.threshold {
                threshold
            } else {
                let calculated = motif.calculate_score(bg, t_offset, renormalize_rows)?;
                if calculated == 0.0 {
                    eprintln!(
                        "Warning: HOMER threshold for motif '{}' was clipped to 0 (t_offset={}); scanning may be very permissive",
                        motif.id, t_offset
                    );
                }
                calculated
            }
        } else {
            let calculated = motif.calculate_score(bg, t_offset, renormalize_rows)?;
            if calculated == 0.0 {
                eprintln!(
                    "Warning: HOMER threshold for motif '{}' was clipped to 0 (t_offset={}); scanning may be very permissive",
                    motif.id, t_offset
                );
            }
            calculated
        };
        writeln!(
            writer,
            ">{}\t{}\t{:.6}\t0\t0\t0",
            motif.id, motif.description, score
        )?;
        for row in &motif.matrix {
            let values = if renormalize_rows {
                renormalized(row)
            } else {
                row.clone()
            };
            let line: Vec<String> = values.iter().map(|v| format!("{:.6}", v)).collect();
            writeln!(writer, "{}", line.join("\t"))?;
        }
    }
    Ok(())
}

pub fn write_meme<W: Write>(
    writer: &mut W,
    motifs: &[Motif],
    nsites: Option<usize>,
    evalue: Option<f64>,
    renormalize_rows: bool,
) -> Result<(), MotifError> {
    let mut header_printed = false;
    let mut header_alphabet = String::new();
    for motif in motifs {
        if !header_printed {
            header_alphabet = motif.alphabet.clone();
            writeln!(writer, "MEME version 4\n")?;
            writeln!(writer, "ALPHABET= {}\n", motif.alphabet)?;
            if motif.alphabet == "ACGT" || motif.alphabet == "ACGU" {
                writeln!(writer, "strands: + -\n")?;
            }
            writeln!(writer, "Background letter frequencies")?;
            writeln!(writer, "{}\n", background_line(&motif.alphabet))?;
            header_printed = true;
        } else if motif.alphabet != header_alphabet {
            eprintln!(
                "Warning: skipping motif '{}' with alphabet {} (header uses {})",
                motif.id, motif.alphabet, header_alphabet
            );
            continue;
        }

        let expected_cols = alphabet_letters(&motif.alphabet).chars().count();
        let width = motif.matrix.len();
        let motif_nsites = nsites.or(motif.nsites).unwrap_or(20);
        let motif_evalue = evalue.or(motif.evalue);
        let evalue_value = motif_evalue
            .map(|v| format!("{:.6}", v))
            .unwrap_or_else(|| "0".to_string());
        writeln!(writer, "MOTIF {} {}", motif.id, motif.description)?;
        writeln!(writer)?;
        writeln!(
            writer,
            "letter-probability matrix: alength= {} w= {} nsites= {} E= {}",
            expected_cols, width, motif_nsites, evalue_value
        )?;
        for row in &motif.matrix {
            let values = if renormalize_rows {
                renormalized(row)
            } else {
                row.clone()
            };
            let line: Vec<String> = values.iter().map(|v| format!("{:.6}", v)).collect();
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
        writeln!(writer, "      \"id\": {},", json_string(&motif.id))?;
        writeln!(
            writer,
            "      \"description\": {},",
            json_string(&motif.description)
        )?;
        if !motif.alphabet.is_empty() {
            writeln!(
                writer,
                "      \"alphabet\": {},",
                json_string(&motif.alphabet)
            )?;
        }
        if let Some(threshold) = motif.threshold {
            writeln!(writer, "      \"threshold\": {:.6},", threshold)?;
        }
        if let Some(nsites) = motif.nsites {
            writeln!(writer, "      \"nsites\": {},", nsites)?;
        }
        if let Some(evalue) = motif.evalue {
            writeln!(writer, "      \"evalue\": {:.6},", evalue)?;
        }
        writeln!(writer, "      \"matrix\": [")?;
        for (ri, row) in motif.matrix.iter().enumerate() {
            let vals: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
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

#[cfg(test)]
mod tests {
    use super::{
        is_logodds, logodds_to_prob, read_homer, read_json, read_meme, write_homer, write_json,
        write_meme, MatrixType,
    };
    use std::io::Cursor;

    #[test]
    fn is_logodds_respects_explicit_matrix_type() {
        assert!(!is_logodds(&[0.2, 0.3, 0.3, 0.2], MatrixType::Probability));
        assert!(is_logodds(&[0.2, 0.3, 0.3, 0.2], MatrixType::Logodds));
        assert!(is_logodds(&[2.0, 0.0, -1.0, -2.0], MatrixType::Auto));
        assert!(!is_logodds(&[0.25, 0.25, 0.25, 0.25], MatrixType::Auto));
    }

    #[test]
    fn logodds_to_prob_normalizes_rows() {
        let row = logodds_to_prob(&[2.0, 0.0, -1.0, -2.0], 0.01, &[0.25]).unwrap();
        let sum: f64 = row.iter().sum();

        assert!((sum - 1.0).abs() < 1e-9);
        assert!(row[0] > row[1]);
        assert!(row[1] > row[2]);
        assert!(row[2] > row[3]);
    }

    #[test]
    fn logodds_to_prob_rejects_background_vector_sum_mismatch() {
        let err =
            logodds_to_prob(&[2.0, 0.0, -1.0, -2.0], 0.01, &[0.50, 0.50, 0.50, 0.50]).unwrap_err();

        assert!(err.to_string().contains("sum to 1.0"));
    }

    #[test]
    fn read_meme_keeps_alphabet_width_when_alength_conflicts() {
        let source = "MEME version 4\n\
\n\
ALPHABET= PROTEIN\n\
\n\
MOTIF P0001 protein_motif\n\
\n\
letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 0\n\
  0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05\n\
  0.10 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.00 0.05\n\
//\n";

        let motifs = read_meme(Cursor::new(source), None).unwrap();

        assert_eq!(motifs.len(), 1);
        assert_eq!(motifs[0].alphabet, "PROTEIN");
        assert_eq!(motifs[0].matrix.len(), 2);
        assert!(motifs[0].matrix.iter().all(|row| row.len() == 20));
        assert_eq!(motifs[0].nsites, Some(20));
        assert_eq!(motifs[0].evalue, Some(0.0));
    }

    #[test]
    fn write_meme_omits_strands_for_protein() {
        let source = format!(
            "ALPHABET= PROTEIN\nMOTIF P1\nletter-probability matrix:\n{}\n",
            vec!["0.05"; 20].join(" ")
        );
        let motifs = read_meme(Cursor::new(source), None).unwrap();
        let mut out = Vec::new();

        write_meme(&mut out, &motifs, None, None, false).unwrap();

        let rendered = String::from_utf8(out).unwrap();
        assert!(rendered.contains("ALPHABET= PROTEIN"));
        assert!(!rendered.contains("strands: + -"));
        assert!(rendered.contains("A 0.05"));
    }

    #[test]
    fn write_meme_overrides_metadata_and_renormalizes_rows() {
        let source = ">M1\tdesc\t1\t0\t0\t0\n0.20\t0.20\t0.20\t0.20\n";
        let motifs = read_homer(
            Cursor::new(source),
            0.01,
            MatrixType::Probability,
            "ACGT",
            &[0.25],
        )
        .unwrap();
        let mut out = Vec::new();

        write_meme(&mut out, &motifs, Some(4000), Some(0.00001), true).unwrap();

        let rendered = String::from_utf8(out).unwrap();
        assert!(rendered.contains("nsites= 4000 E= 0.000010"));
        assert!(rendered.contains("0.250000  0.250000  0.250000  0.250000"));
    }

    #[test]
    fn write_meme_uses_motif_metadata_without_cli_override() {
        let source = "MOTIF M1 desc\n\
letter-probability matrix: alength= 4 w= 1 nsites= 123 E= 0.5\n\
0.25 0.25 0.25 0.25\n";
        let motifs = read_meme(Cursor::new(source), None).unwrap();
        let mut out = Vec::new();

        write_meme(&mut out, &motifs, None, None, false).unwrap();

        let rendered = String::from_utf8(out).unwrap();
        assert!(rendered.contains("nsites= 123 E= 0.500000"));
    }

    #[test]
    fn json_round_trip_preserves_threshold_nsites_and_evalue_metadata() {
        let source = "MOTIF M1 desc\n\
letter-probability matrix: alength= 4 w= 1 nsites= 123 E= 0.5\n\
0.25 0.25 0.25 0.25\n";
        let mut motifs = read_meme(Cursor::new(source), None).unwrap();
        motifs[0].threshold = Some(7.25);
        let mut json = Vec::new();

        write_json(&mut json, &motifs).unwrap();
        let rendered = String::from_utf8(json).unwrap();
        let parsed = read_json(
            Cursor::new(rendered.as_bytes()),
            0.01,
            MatrixType::Probability,
            &[0.25],
        )
        .unwrap();

        assert!(rendered.contains("\"threshold\": 7.250000"));
        assert!(rendered.contains("\"nsites\": 123"));
        assert!(rendered.contains("\"evalue\": 0.500000"));
        assert_eq!(parsed[0].threshold, Some(7.25));
        assert_eq!(parsed[0].nsites, Some(123));
        assert_eq!(parsed[0].evalue, Some(0.5));
    }

    #[test]
    fn write_homer_can_keep_source_threshold() {
        let source = ">KEEP\tdesc\t9.5\t0\t0\t0\n0.25\t0.25\t0.25\t0.25\n";
        let motifs = read_homer(
            Cursor::new(source),
            0.01,
            MatrixType::Probability,
            "ACGT",
            &[0.25],
        )
        .unwrap();
        let mut out = Vec::new();

        write_homer(&mut out, &motifs, &[0.25], 4.0, true, false).unwrap();

        let rendered = String::from_utf8(out).unwrap();
        assert!(rendered.starts_with(">KEEP\tdesc\t9.500000\t"));
    }
}
