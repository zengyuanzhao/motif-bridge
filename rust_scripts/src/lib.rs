pub mod error;
pub mod io;

pub use error::MotifError;
pub use io::MatrixType;
use std::f64;

#[derive(Clone, Debug)]
pub struct Motif {
    pub id: String,
    pub description: String,
    pub matrix: Vec<Vec<f64>>,
    pub alphabet: String,
    pub threshold: Option<f64>,
    pub nsites: Option<usize>,
    pub evalue: Option<f64>,
}

impl Motif {
    pub fn new(id: String, description: String, matrix: Vec<Vec<f64>>, alphabet: String) -> Self {
        Self {
            id,
            description,
            matrix,
            alphabet,
            threshold: None,
            nsites: None,
            evalue: None,
        }
    }

    pub fn with_threshold(
        id: String,
        description: String,
        matrix: Vec<Vec<f64>>,
        alphabet: String,
        threshold: Option<f64>,
    ) -> Self {
        Self::with_metadata(id, description, matrix, alphabet, threshold, None, None)
    }

    pub fn with_metadata(
        id: String,
        description: String,
        matrix: Vec<Vec<f64>>,
        alphabet: String,
        threshold: Option<f64>,
        nsites: Option<usize>,
        evalue: Option<f64>,
    ) -> Self {
        Self {
            id,
            description,
            matrix,
            alphabet,
            threshold,
            nsites,
            evalue,
        }
    }

    pub fn calculate_score(
        &self,
        bg: &[f64],
        t_offset: f64,
        renormalize: bool,
    ) -> Result<f64, MotifError> {
        let mut raw = 0.0;
        for row in &self.matrix {
            if row.is_empty() {
                continue;
            }
            let values = if renormalize {
                renormalized(row)
            } else {
                row.clone()
            };
            let backgrounds = background_values(bg, values.len())?;
            let mut best_idx = 0;
            for (i, &value) in values.iter().enumerate() {
                if value > values[best_idx] {
                    best_idx = i;
                }
            }
            let max_p = values[best_idx];
            if max_p > 0.0 {
                raw += (max_p / backgrounds[best_idx]).log2();
            }
        }
        Ok((raw - t_offset).max(0.0))
    }

    pub fn calculate_ic(&self) -> Vec<f64> {
        let is_protein = self.alphabet == "PROTEIN";
        let max_ic = if is_protein { 20.0_f64.log2() } else { 2.0 };

        self.matrix
            .iter()
            .map(|row| {
                let mut h = 0.0;
                for &p in row {
                    if p > 0.0 {
                        h -= p * p.log2();
                    }
                }
                let ic = max_ic - h;
                if ic < 0.0 {
                    0.0
                } else {
                    ic
                }
            })
            .collect()
    }

    pub fn total_ic(&self) -> f64 {
        self.calculate_ic().iter().sum()
    }

    pub fn trim_edges(&mut self, threshold: f64) {
        let ic_list = self.calculate_ic();
        let mut start = 0;
        while start < ic_list.len() && ic_list[start] < threshold {
            start += 1;
        }
        let mut end = ic_list.len();
        while end > start && ic_list[end - 1] < threshold {
            end -= 1;
        }
        if start < end {
            self.matrix = self.matrix[start..end].to_vec();
        } else {
            self.matrix.clear();
        }
    }

    pub fn reverse_complement(&mut self) -> Result<(), String> {
        if self.alphabet != "ACGT" && self.alphabet != "ACGU" {
            return Err(format!(
                "Reverse complement not supported for alphabet: {}",
                self.alphabet
            ));
        }
        let mut rc_matrix = Vec::with_capacity(self.matrix.len());
        for row in self.matrix.iter().rev() {
            if row.len() == 4 {
                rc_matrix.push(vec![row[3], row[2], row[1], row[0]]);
            } else {
                return Err("Row length is not 4".to_string());
            }
        }
        self.matrix = rc_matrix;
        self.id.push_str("_RC");
        Ok(())
    }

    // Direct-print helpers removed; use io::write_homer/write_meme for output.
}

pub fn background_values(bg: &[f64], width: usize) -> Result<Vec<f64>, MotifError> {
    if bg.len() == 1 {
        if bg[0] <= 0.0 || bg[0] > 1.0 {
            return Err(MotifError::Operation(
                "background values must be in (0, 1]".to_string(),
            ));
        }
        return Ok(vec![bg[0]; width]);
    }
    if bg.len() != width {
        return Err(MotifError::Operation(format!(
            "background length {} does not match row width {}",
            bg.len(),
            width
        )));
    }
    if bg.iter().any(|v| *v <= 0.0 || *v > 1.0) {
        return Err(MotifError::Operation(
            "background values must be in (0, 1]".to_string(),
        ));
    }
    let total: f64 = bg.iter().sum();
    if (total - 1.0).abs() > 1e-3 {
        return Err(MotifError::Operation(format!(
            "background vector must sum to 1.0, got {:.6}",
            total
        )));
    }
    Ok(bg.to_vec())
}

pub fn renormalized(row: &[f64]) -> Vec<f64> {
    let total: f64 = row.iter().sum();
    if total > 0.0 {
        row.iter().map(|v| v / total).collect()
    } else {
        row.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::Motif;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-9,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn calculate_ic_dna_conserved_and_uniform_rows() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.25, 0.25, 0.25, 0.25]],
            "ACGT".to_string(),
        );

        let ic = motif.calculate_ic();

        assert_close(ic[0], 2.0);
        assert_close(ic[1], 0.0);
        assert_close(motif.total_ic(), 2.0);
    }

    #[test]
    fn calculate_ic_protein_uses_twenty_letter_maximum() {
        let mut conserved = vec![0.0; 20];
        conserved[0] = 1.0;
        let uniform = vec![0.05; 20];
        let motif = Motif::new(
            "p1".to_string(),
            "protein".to_string(),
            vec![conserved, uniform],
            "PROTEIN".to_string(),
        );

        let ic = motif.calculate_ic();

        assert_close(ic[0], 20.0_f64.log2());
        assert_close(ic[1], 0.0);
    }

    #[test]
    fn reverse_complement_reverses_rows_and_swaps_columns() {
        let mut motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.1, 0.2, 0.3, 0.4], vec![0.5, 0.6, 0.7, 0.8]],
            "ACGT".to_string(),
        );

        motif.reverse_complement().unwrap();

        assert_eq!(motif.id, "m1_RC");
        assert_eq!(
            motif.matrix,
            vec![vec![0.8, 0.7, 0.6, 0.5], vec![0.4, 0.3, 0.2, 0.1]]
        );
    }

    #[test]
    fn reverse_complement_rejects_protein() {
        let mut motif = Motif::new(
            "p1".to_string(),
            "protein".to_string(),
            vec![vec![1.0; 20]],
            "PROTEIN".to_string(),
        );

        let err = motif.reverse_complement().unwrap_err();

        assert!(err.contains("Reverse complement not supported"));
    }

    #[test]
    fn trim_edges_removes_low_ic_flanks() {
        let mut motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![
                vec![0.25, 0.25, 0.25, 0.25],
                vec![1.0, 0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0, 0.0],
                vec![0.25, 0.25, 0.25, 0.25],
            ],
            "ACGT".to_string(),
        );

        motif.trim_edges(0.5);

        assert_eq!(
            motif.matrix,
            vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn calculate_score_uses_column_specific_background() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.1, 0.2, 0.3, 0.4], vec![0.6, 0.2, 0.1, 0.1]],
            "ACGT".to_string(),
        );

        let score = motif
            .calculate_score(&[0.30, 0.20, 0.20, 0.30], 0.0, false)
            .unwrap();

        assert_close(score, (0.4_f64 / 0.30).log2() + (0.6_f64 / 0.30).log2());
    }

    #[test]
    fn calculate_score_uses_first_column_when_maximum_ties() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.4, 0.1, 0.1, 0.4]],
            "ACGT".to_string(),
        );

        let score = motif
            .calculate_score(&[0.10, 0.20, 0.30, 0.40], 0.0, false)
            .unwrap();

        assert_close(score, (0.4_f64 / 0.10).log2());
    }

    #[test]
    fn calculate_score_rejects_background_width_mismatch() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.25, 0.25, 0.25, 0.25]],
            "ACGT".to_string(),
        );

        let err = motif
            .calculate_score(&[0.25, 0.25], 0.0, false)
            .unwrap_err();

        assert!(err.to_string().contains("background length"));
    }

    #[test]
    fn calculate_score_rejects_invalid_background_values() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.25, 0.25, 0.25, 0.25]],
            "ACGT".to_string(),
        );

        let err = motif
            .calculate_score(&[0.25, 0.25, 0.25, 0.0], 0.0, false)
            .unwrap_err();

        assert!(err.to_string().contains("background values"));
    }

    #[test]
    fn calculate_score_rejects_background_vector_sum_mismatch() {
        let motif = Motif::new(
            "m1".to_string(),
            "dna".to_string(),
            vec![vec![0.25, 0.25, 0.25, 0.25]],
            "ACGT".to_string(),
        );

        let err = motif
            .calculate_score(&[0.50, 0.50, 0.50, 0.50], 0.0, false)
            .unwrap_err();

        assert!(err.to_string().contains("sum to 1.0"));
    }
}
