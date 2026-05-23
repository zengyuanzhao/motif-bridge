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
}

impl Motif {
    pub fn new(id: String, description: String, matrix: Vec<Vec<f64>>, alphabet: String) -> Self {
        Self {
            id,
            description,
            matrix,
            alphabet,
        }
    }

    pub fn calculate_score(&self, bg: f64, t_offset: f64) -> f64 {
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

    pub fn print_homer(&self, bg: f64, t_offset: f64) {
        let score = self.calculate_score(bg, t_offset);
        println!(">{}\t{}\t{:.6}\t0\t0\t0", self.id, self.description, score);
        for row in &self.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            println!("{}", line.join("\t"));
        }
    }

    pub fn print_meme_motif(&self) {
        let expected_cols = if self.alphabet == "PROTEIN" { 20 } else { 4 };
        let width = self.matrix.len();
        println!("MOTIF {} {}", self.id, self.description);
        println!();
        println!(
            "letter-probability matrix: alength= {} w= {} nsites= 20 E= 0",
            expected_cols, width
        );
        for row in &self.matrix {
            let line: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            println!("  {}", line.join("  "));
        }
        println!();
    }
}
