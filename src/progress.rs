use std::io::{self, Write};

/// Progress bar structure
pub struct ProgressBar {
    total: usize,
    current: usize,
    width: usize,
    last_percentage: usize,
}

impl ProgressBar {
    /// Create a new progress bar
    pub fn new(total: usize) -> Self {
        Self {
            total,
            current: 0,
            width: 50, // Progress bar width
            last_percentage: 0,
        }
    }
    
    /// Set progress bar width
    pub fn with_width(mut self, width: usize) -> Self {
        self.width = width;
        self
    }
    
    /// Update progress
    pub fn update(&mut self, current: usize) -> io::Result<()> {
        self.current = current;
        let percentage = if self.total > 0 {
            (current * 100) / self.total
        } else {
            0
        };
        
        // Only update display when percentage changes
        if percentage != self.last_percentage {
            self.display()?;
            self.last_percentage = percentage;
        }
        
        Ok(())
    }
    
    /// Display progress bar
    fn display(&self) -> io::Result<()> {
        let percentage = if self.total > 0 {
            (self.current * 100) / self.total
        } else {
            0
        };
        
        let filled_width = if self.total > 0 {
            (self.current * self.width) / self.total
        } else {
            0
        };
        
        let bar = "█".repeat(filled_width);
        let empty = "░".repeat(self.width - filled_width);
        
        // Clear current line and display progress bar
        print!("\r[{}] {}% ({}/{})", 
               bar + &empty, 
               percentage, 
               self.current, 
               self.total);
        io::stdout().flush()?;
        
        Ok(())
    }
    
    /// Finish progress bar
    pub fn finish(&mut self) -> io::Result<()> {
        self.current = self.total;
        self.display()?;
        println!(); // New line
        Ok(())
    }
    
    /// Display progress bar with description
    pub fn with_description(mut self, description: &str) -> io::Result<()> {
        let percentage = if self.total > 0 {
            (self.current * 100) / self.total
        } else {
            0
        };
        
        let filled_width = if self.total > 0 {
            (self.current * self.width) / self.total
        } else {
            0
        };
        
        let bar = "█".repeat(filled_width);
        let empty = "░".repeat(self.width - filled_width);
        
        // Clear current line and display progress bar with description
        print!("\r{} [{}] {}% ({}/{})", 
               description,
               bar + &empty, 
               percentage, 
               self.current, 
               self.total);
        io::stdout().flush()?;
        
        Ok(())
    }
}

/// Simple progress displayer
pub struct SimpleProgress {
    total: usize,
    current: usize,
    last_percentage: usize,
}

impl SimpleProgress {
    /// Create a new simple progress displayer
    pub fn new(total: usize) -> Self {
        Self {
            total,
            current: 0,
            last_percentage: 0,
        }
    }
    
    /// Update progress (refresh on each call to avoid staying at fixed count for long time)
    pub fn update(&mut self, current: usize) -> io::Result<()> {
        self.current = current;
        let percentage = if self.total > 0 {
            (current * 100) / self.total
        } else {
            0
        };

        print!("\r[Progressing] {}/{} ({}%)", 
               self.current, 
               self.total, 
               percentage);
        io::stdout().flush()?;
        self.last_percentage = percentage;
        
        Ok(())
    }
    
    /// Finish progress display
    pub fn finish(&mut self) -> io::Result<()> {
        self.current = self.total;
        println!("\r[Progressing] {}/{} (100%)", self.total, self.total);
        io::stdout().flush()?;
        Ok(())
    }
}

/// Descriptive progress displayer
pub struct DescriptiveProgress {
    total: usize,
    current: usize,
    description: String,
    last_percentage: usize,
}

impl DescriptiveProgress {
    /// Create a new descriptive progress displayer
    pub fn new(total: usize, description: &str) -> Self {
        Self {
            total,
            current: 0,
            description: description.to_string(),
            last_percentage: 0,
        }
    }
    
    /// Update progress
    pub fn update(&mut self, current: usize) -> io::Result<()> {
        self.current = current;
        let percentage = if self.total > 0 {
            (current * 100) / self.total
        } else {
            0
        };
        
        // Only update display when percentage changes
        if percentage != self.last_percentage {
            print!("\r{}: {}/{} ({}%)", 
                   self.description,
                   self.current, 
                   self.total, 
                   percentage);
            io::stdout().flush()?;
            self.last_percentage = percentage;
        }
        
        Ok(())
    }
    
    /// Finish progress display
    pub fn finish(&mut self) -> io::Result<()> {
        self.current = self.total;
        print!("\r{}: {}/{} (100%)", 
               self.description,
               self.total, 
               self.total);
        println!(); // New line
        Ok(())
    }
} 

/// Format time as "xx h xx m xx.xxx s" format
pub fn format_time_used(elapsed: std::time::Duration) -> String {
    let total_secs = elapsed.as_secs_f64();
    let hours = (total_secs / 3600.0) as u64;
    let minutes = ((total_secs % 3600.0) / 60.0) as u64;
    let seconds = total_secs % 60.0;
    
    if hours > 0 {
        format!("[Time used] {:02} h {:02} m {:05.3} s", hours, minutes, seconds)
    } else if minutes > 0 {
        format!("[Time used] {:02} m {:05.3} s", minutes, seconds)
    } else {
        format!("[Time used] {:05.3} s", seconds)
    }
}