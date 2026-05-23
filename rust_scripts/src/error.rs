use std::fmt;
use std::io;

#[derive(Debug)]
pub enum MotifError {
    Io(io::Error),
    Parse(String),
    Operation(String),
}

impl fmt::Display for MotifError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MotifError::Io(err) => write!(f, "IO error: {}", err),
            MotifError::Parse(msg) => write!(f, "Parse error: {}", msg),
            MotifError::Operation(msg) => write!(f, "Operation error: {}", msg),
        }
    }
}

impl std::error::Error for MotifError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MotifError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for MotifError {
    fn from(err: io::Error) -> Self {
        MotifError::Io(err)
    }
}

impl From<serde_json::Error> for MotifError {
    fn from(err: serde_json::Error) -> Self {
        MotifError::Parse(format!("JSON parsing error: {}", err))
    }
}
