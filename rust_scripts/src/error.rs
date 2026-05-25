use std::io;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MotifError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Operation error: {0}")]
    Operation(String),
    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),
}
