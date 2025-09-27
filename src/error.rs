use thiserror::Error;
#[derive(Error, Debug)]
pub enum AppError {
    #[error("Failed to read BAM file: {0}")]
    BamRead(#[from] rust_htslib::errors::Error),
    #[error("Failed to parse configuration: {0}")]
    ConfigParse(#[from] config::ConfigError),
}