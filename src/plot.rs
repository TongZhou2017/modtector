use chrono;
use plotters::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/// Parse f64 value from string, supporting "NaN" string
fn parse_f64_allow_nan(s: &str) -> Option<f64> {
    let trimmed = s.trim();
    if trimmed.eq_ignore_ascii_case("nan") {
        Some(f64::NAN)
    } else {
        trimmed.parse::<f64>().ok()
    }
}

/// SVG drawing configuration structure
#[derive(Debug, Clone)]
pub struct SvgDrawConfig {
    /// Color circle radius
    pub circle_radius: f64,
    /// Color circle stroke width
    pub circle_stroke_width: f64,
    /// Whether to draw stroke
    pub circle_has_stroke: bool,
    /// Circle type: true means filled circle, false means hollow ring
    pub circle_filled: bool,
    /// Font color setting: Some(color) means use specified color, None means use template default color
    pub font_color: Option<String>,
    /// Layer order: true means color circles before text, false means color circles after text
    pub circles_before_text: bool,
    /// Whether to draw color bars
    pub draw_color_bars: bool,
    /// Legend each color segment width (when horizontal layout)
    pub legend_item_width: f64,
    /// Legend each color segment height (when vertical layout)
    pub legend_item_height: f64,
}

impl Default for SvgDrawConfig {
    fn default() -> Self {
        Self {
            circle_radius: 2.8,
            circle_stroke_width: 0.2,
            circle_has_stroke: false,
            circle_filled: false,       // Default hollow ring
            font_color: None,           // Default use template color
            circles_before_text: false, // Default color circles after text, text on top layer
            draw_color_bars: true,
            legend_item_width: 30.0,    // Default legend each color segment width
            legend_item_height: 15.0,   // Default legend each color segment height
        }
    }
}

/// Color circle information
#[derive(Debug, Clone)]
pub struct ColorCircle {
    pub x: f64,
    pub y: f64,
    pub color: String,
    pub radius: f64,
    pub stroke_width: f64,
    pub has_stroke: bool,
    pub position: Option<u32>, // Position for interactive visualization
    pub base: Option<char>,     // Base type for interactive visualization
    pub reactivity: Option<f64>, // Reactivity value for interactive visualization
}

/// SVG element information
#[derive(Debug, Clone)]
pub struct SvgElement {
    pub element_type: String, // "circle", "text", "line", etc.
    pub content: String,
    pub layer: u32, // Layer order, smaller numbers are at lower layers
}

/// Parse SVG template to extract position information where color circles need to be added
pub fn parse_svg_template_for_circles(
    svg_content: &str,
    pos_to_reactivity: &HashMap<u32, (char, f64)>,
    color_ranges: &BaseColorRanges,
    config: &SvgDrawConfig,
    base_shifts: &HashMap<char, i32>,
    template_pos2base: &HashMap<u32, char>,
) -> Result<Vec<ColorCircle>, Box<dyn Error>> {
    let position_regex = Regex::new(r"<title>(\d+)")?;
    // Priority: extract position from "position.label in template: XXX" if available
    let template_position_regex = Regex::new(r"position\.label in template: (\d+)")?;
    let x_regex = Regex::new(r#"x="([^"]+)""#)?;
    let y_regex = Regex::new(r#"y="([^"]+)""#)?;

    let mut circles = Vec::new();

    // Extract base from SVG template line to determine shift
    let base_regex = Regex::new(r"<text[^>]*>([ATCGU])</text>")?;
    
    for template_line in svg_content.lines() {
        if template_line.contains("title>") && !template_line.contains("numbering-label") {
            // Priority: use template position if available, otherwise use title position
            let svg_position_opt = if let Some(template_caps) = template_position_regex.captures(template_line) {
                template_caps[1].parse::<u32>().ok()
            } else if let Some(title_caps) = position_regex.captures(template_line) {
                title_caps[1].parse::<u32>().ok()
            } else {
                None
            };
            
            if let Some(svg_position) = svg_position_opt {
                // CRITICAL FIX: Use template_pos2base mapping (from "position.label in template: POS.BASE")
                // This is the most reliable source of truth for base types in the SVG template
                // Do NOT use "nearby search <text>" method, as it can match wrong elements
                let base_from_template = template_pos2base.get(&svg_position).copied();
                
                // Calculate CSV position with shift
                // Shift formula: CSV_position = SVG_position + shift
                // But if SVG_position is already a genomic coordinate (>= 100000), use it directly
                let csv_position = if let Some(base) = base_from_template {
                    let normalized_base = normalize_base(base);
                    // Check if SVG position is already a genomic coordinate (>= 100000)
                    // If so, use it directly without shift
                    if svg_position >= 100000 {
                        svg_position // Already genomic coordinate, use directly
                    } else if let Some(&shift) = base_shifts.get(&normalized_base) {
                        // Relative position, apply shift
                        let shifted_pos = (svg_position as i32) + shift;
                        if shifted_pos > 0 {
                            shifted_pos as u32
                        } else {
                            svg_position // Fallback to original position
                        }
                    } else {
                        svg_position // No shift for this base
                    }
                } else {
                    // If no base found, check if position is genomic coordinate
                    if svg_position >= 100000 {
                        svg_position // Already genomic coordinate, use directly
                    } else {
                        svg_position // No base found, use original position
                    }
                };
                
                // CRITICAL: Only process if we have a valid base from template
                // If template doesn't have this position, skip it to avoid wrong matches
                if let Some(svg_base) = base_from_template {
                    let normalized_svg_base = normalize_base(svg_base);
                    
                    if csv_position > 0 {
                        if let Some((csv_base, reactivity)) = pos_to_reactivity.get(&csv_position) {
                            // Normalize bases for comparison (T and U are equivalent)
                            let normalized_csv_base = normalize_base(*csv_base);
                            
                            // Check if bases match (considering T/U equivalence)
                            if normalized_csv_base != normalized_svg_base {
                                continue; // Skip if bases don't match
                            }
                            
                            // CRITICAL FIX: Use SVG template's base type (normalized_svg_base) for data-base attribute
                            // This ensures consistency with what's actually displayed in the SVG
                            // The SVG template's base is the source of truth for visualization
                            let base_for_visualization = normalized_svg_base;
                            
                            let color_ranges_for_base = match normalized_csv_base {
                                'A' => &color_ranges.a,
                                'T' => &color_ranges.t,
                                'C' => &color_ranges.c,
                                'G' => &color_ranges.g,
                                _ => continue,
                            };
                            let color = determine_color(*reactivity, color_ranges_for_base);

                            // Extract x and y coordinates
                            if let (Some(x_caps), Some(y_caps)) = (
                                x_regex.captures(template_line),
                                y_regex.captures(template_line),
                            ) {
                                if let (Ok(x), Ok(y)) =
                                    (x_caps[1].parse::<f64>(), y_caps[1].parse::<f64>())
                                {
                                    circles.push(ColorCircle {
                                        x,
                                        y,
                                        color,
                                        radius: config.circle_radius,
                                        stroke_width: config.circle_stroke_width,
                                        has_stroke: config.circle_has_stroke
                                            || !config.circle_filled, // Force stroke when hollow ring
                                        position: Some(csv_position),
                                        base: Some(base_for_visualization), // CRITICAL: Use SVG template's base, not CSV base
                                        reactivity: Some(*reactivity), // Store reactivity value
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(circles)
}

/// Draw color circles to SVG
pub fn draw_color_circles(
    writer: &mut BufWriter<File>,
    circles: &[ColorCircle],
    config: &SvgDrawConfig,
) -> Result<(), Box<dyn Error>> {
    for circle in circles {
        // Use inline style to ensure priority over template CSS (e.g., circle { stroke: black })
        let style_string = if config.circle_filled {
            if circle.has_stroke {
                format!(
                    "fill: {}; stroke: {}; stroke-width: {};",
                    circle.color, circle.color, circle.stroke_width
                )
            } else {
                format!("fill: {}; stroke: none;", circle.color)
            }
        } else {
            // Hollow ring: transparent fill + colored stroke
            format!(
                "fill: none; stroke: {}; stroke-width: {};",
                circle.color, circle.stroke_width
            )
        };

        // Add title element and data attributes for interactive visualization if position is available
        // CRITICAL: Use unique physical ID for atomic DOM access (eliminates selector ambiguity)
        if let Some(position) = circle.position {
            let base_attr = if let Some(base) = circle.base {
                format!(r#" data-base="{}""#, base)
            } else {
                String::new()
            };
            let reactivity_attr = if let Some(reactivity) = circle.reactivity {
                format!(r#" data-reactivity="{:.6}""#, reactivity)
            } else {
                String::new()
            };
            // Generate unique physical ID: node_{position} for direct DOM access
            let node_id = format!("node_{}", position);
            writeln!(
                writer,
                r#"<g><title>{}</title><circle id="{}" class="rna-base-circle" cx="{}" cy="{}" r="{}" style="{}" data-position="{}"{}{} /></g>"#,
                position, node_id, circle.x, circle.y, circle.radius, style_string, position, base_attr, reactivity_attr
            )?;
        } else {
        writeln!(
            writer,
            r#"<g><circle cx="{}" cy="{}" r="{}" style="{}" /></g>"#,
            circle.x, circle.y, circle.radius, style_string
        )?;
        }
    }
    Ok(())
}

/// Draw SVG skeleton and text
pub fn draw_svg_skeleton_and_text(
    writer: &mut BufWriter<File>,
    svg_content: &str,
) -> Result<(), Box<dyn Error>> {
    for template_line in svg_content.lines() {
        if template_line.contains("</svg>") {
            break; // Stop before closing tag
        }
        writeln!(writer, "{}", template_line)?;
    }
    Ok(())
}

/// Draw SVG skeleton (without text)
pub fn draw_svg_skeleton_without_text(
    writer: &mut BufWriter<File>,
    svg_content: &str,
) -> Result<(), Box<dyn Error>> {
    let mut skipped_lines = 0;
    let mut total_lines = 0;

    for template_line in svg_content.lines() {
        if template_line.contains("</svg>") {
            break; // Stop before closing tag
        }
        total_lines += 1;

        // Skip lines containing text elements, but keep other content
        if template_line.contains("<text") && template_line.contains("</text>") {
            // This is a complete text element line, skip it
            skipped_lines += 1;
            // DEBUG output removed to reduce noise
            continue;
        }
        writeln!(writer, "{}", template_line)?;
    }

    // Summary output removed to reduce noise (uncomment for debugging if needed)
    // eprintln!(
    //     "draw_svg_skeleton_without_text: processed {} lines, skipped {} text lines",
    //     total_lines, skipped_lines
    // );
    Ok(())
}

/// Only draw text part in SVG
pub fn draw_svg_text_only(
    writer: &mut BufWriter<File>,
    svg_content: &str,
    config: &SvgDrawConfig,
) -> Result<(), Box<dyn Error>> {
    let mut written_lines = 0;
    for template_line in svg_content.lines() {
        if template_line.contains("</svg>") {
            break; // Stop before closing tag
        }
        
        // Only write lines containing text elements
        if template_line.contains("<text") {
        let processed_line = if let Some(ref font_color) = config.font_color {
            replace_text_color_classes(template_line, font_color)
        } else {
            template_line.to_string()
        };
        writeln!(writer, "{}", processed_line)?;
            written_lines += 1;
    }
    }
    // Summary output removed to reduce noise (uncomment for debugging if needed)
    // eprintln!("draw_svg_text_only: wrote {} text lines", written_lines);
    Ok(())
}

/// Replace color classes in text lines
fn replace_text_color_classes(line: &str, font_color: &str) -> String {
    // If it's a <text ...> element line: remove class="..." and force fill="..."
    if line.contains("<text") {
        let class_re = Regex::new(r#"\sclass=\"[^\"]*\""#).unwrap();
        let mut s = class_re.replace_all(line, "").to_string();

        let fill_re = Regex::new(r#"\sfill=\"[^\"]*\""#).unwrap();
        if fill_re.is_match(&s) {
            s = fill_re
                .replace_all(&s, format!(r#" fill=\"{}\""#, font_color))
                .to_string();
        } else {
            if let Some(pos) = s.find('>') {
                s = format!(
                    "{} {}{}",
                    &s[..pos],
                    &format!(r#"fill=\"{}\""#, font_color),
                    &s[pos..]
                );
            }
        }
        return s;
    }

    // If it's a style line text { fill: ... }, update the fill value
    if line.contains("text {") && line.contains("fill:") {
        let re = Regex::new(r"text\s*\{\s*fill:\s*[^;]+;").unwrap();
        return re
            .replace(line, &format!("text {{fill: {};", font_color))
            .to_string();
    }

    // Other lines unchanged
    line.to_string()
}

/// Draw color bars
pub fn draw_color_bars(
    writer: &mut BufWriter<File>,
    base_groups: &HashMap<char, Vec<SvgReactivityData>>,
    color_ranges: &BaseColorRanges,
    bases_filter: &str,
    config: &SvgDrawConfig,
    svg_width: Option<f64>,
    svg_height: Option<f64>,
) -> Result<(), Box<dyn Error>> {
    // Count total number of bases to calculate total legend height
    let total_bases = bases_filter.chars()
        .filter(|&base| base_groups.contains_key(&base))
        .count();
    
    for (base_idx, base) in bases_filter.chars().enumerate() {
        if let Some(_data_list) = base_groups.get(&base) {
            let color_ranges_for_base = match base {
                'A' => &color_ranges.a,
                'T' => &color_ranges.t,
                'C' => &color_ranges.c,
                'G' => &color_ranges.g,
                _ => continue,
            };
            draw_color_bar(writer, color_ranges_for_base, base, base_idx, config, svg_width, svg_height, total_bases)?;
        }
    }
    Ok(())
}

/// Main SVG drawing function, supports custom layer order
pub fn draw_svg_with_custom_layers(
    reactivity_file: &str,
    svg_template: &str,
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    color_ranges: Option<BaseColorRanges>,
    config: SvgDrawConfig,
    ref_sequence_file: Option<&str>,
    max_shift: u32,
) -> Result<(), Box<dyn Error>> {
    // DEBUG output removed to reduce noise
    // eprintln!("=== draw_svg_with_custom_layers called ===");
    // eprintln!("Using config: {:?}", config);
    // eprintln!(
    //     "DEBUG: draw_svg_with_custom_layers called with circles_before_text = {}",
    //     config.circles_before_text
    // );

    // Write debug info to a file
    if let Ok(mut f) = std::fs::File::create("/tmp/debug_draw_svg.txt") {
        use std::io::Write;
        writeln!(
            f,
            "draw_svg_with_custom_layers called with circles_before_text = {}",
            config.circles_before_text
        )
        .ok();
    }

    let reactivity_data = load_svg_reactivity_data(reactivity_file, signal_type, strand_filter)?;
    // DEBUG output removed to reduce noise
    // eprintln!(
    //     "DEBUG: About to process {} reactivity data points",
    //     reactivity_data.len()
    // );
    let svg_content = std::fs::read_to_string(svg_template)?;
    let color_ranges = color_ranges.unwrap_or_default();

    // Filter data by bases and group by base type
    let bases: Vec<char> = bases_filter.chars().collect();
    let mut base_groups: HashMap<char, Vec<SvgReactivityData>> = HashMap::new();
    for data in &reactivity_data {
        if bases.contains(&data.base) {
            base_groups
                .entry(data.base)
                .or_insert_with(Vec::new)
                .push(data.clone());
        }
    }

    // Create position to reactivity mapping for all bases
    // When multiple signals exist for the same position, prefer non-NaN values
    let mut pos_to_reactivity: HashMap<u32, (char, f64)> = HashMap::new();
    for data in &reactivity_data {
        if bases.contains(&data.base) {
            let existing = pos_to_reactivity.get(&data.position);
            let should_insert = match existing {
                None => true, // No existing entry, insert
                Some((_, existing_reactivity)) => {
                    // Prefer non-NaN values: if existing is NaN and new is not, replace
                    // If existing is valid and new is NaN, keep existing
                    // If both are valid, replace (allows mutation to override stop)
                    if existing_reactivity.is_nan() && !data.reactivity.is_nan() {
                        true // Replace NaN with valid value
                    } else if !existing_reactivity.is_nan() && data.reactivity.is_nan() {
                        false // Keep existing valid value, don't replace with NaN
                    } else {
                        // Both are valid or both are NaN: replace (allows later signals to override)
                        true
                    }
                }
            };
            if should_insert {
                pos_to_reactivity.insert(data.position, (data.base, data.reactivity));
            }
        }
    }

    // Load reference sequence and calculate shifts if provided
    // Shift calculation: CSV_position = SVG_position + shift
    // We need to find the shift that maps SVG RNA positions to CSV genomic positions
    let mut base_shifts: HashMap<char, i32> = HashMap::new();
    if let Some(ref_file) = ref_sequence_file {
        // DEBUG output removed to reduce noise
        // eprintln!("Loading reference sequence from: {}", ref_file);
        let ref_pos2base = load_reference_sequence(ref_file)?;
        // DEBUG output removed to reduce noise
        // eprintln!("Loaded {} positions from reference sequence", ref_pos2base.len());
        
        // Extract SVG positions and bases from template
        let svg_pos2base = extract_svg_positions_and_bases(&svg_content)?;
        // DEBUG output removed to reduce noise
        // eprintln!("Extracted {} positions from SVG template", svg_pos2base.len());
        
                // Calculate best shift for each base type
                // Shift maps SVG RNA position to CSV genomic position: CSV_pos = SVG_pos + shift
                // Only use relative positions (< 100000) for shift calculation, as positions >= 100000 are already genomic coordinates
                for base in &bases {
                    // Filter data by normalized base (T and U are equivalent)
                    let base_data: Vec<&SvgReactivityData> = reactivity_data
                        .iter()
                        .filter(|d| normalize_base(d.base) == *base)
                        .collect();
                    
                    if !base_data.is_empty() {
                        // Find shift that maximizes matches between SVG positions and CSV positions
                        // Only consider relative positions (< 100000) for shift calculation
                        let min_csv_pos = base_data.iter().map(|d| d.position).min().unwrap_or(0);
                        // Get relative SVG positions for THIS base type only
                        let relative_svg_positions_for_base: Vec<u32> = svg_pos2base.iter()
                            .filter(|(pos, svg_base)| {
                                **pos < 100000 && normalize_base(**svg_base) == *base
                            })
                            .map(|(pos, _)| *pos)
                            .collect();
                        
                        let min_svg_pos = relative_svg_positions_for_base.iter().min().copied().unwrap_or(0);
                        let estimated_shift = min_csv_pos as i32 - min_svg_pos as i32;
                        
                        eprintln!("Base {}: min_csv_pos={}, min_svg_pos (relative)={}, estimated_shift={}", 
                                 base, min_csv_pos, min_svg_pos, estimated_shift);
                        
                        // Search around the estimated shift
                        let search_range = max_shift.max(1000) as i32; // Use larger range for large shifts
                        let mut best_shift = estimated_shift;
                        let mut best_match_count = 0;
                        let mut best_total_checked = 0;
                        
                        for shift_offset in -search_range..=search_range {
                            let shift = estimated_shift + shift_offset;
                            let mut match_count = 0;
                            let mut total_checked = 0;
                            
                            // Check matches: for each SVG position, check if CSV position (SVG_pos + shift) has matching base
                            // Only check relative positions (< 100000) for shift calculation
                            // Positions >= 100000 are already genomic coordinates and will be used directly
                            for (svg_pos, svg_base) in &svg_pos2base {
                                let normalized_svg_base = normalize_base(*svg_base);
                                if normalized_svg_base == *base {
                                    // Skip genomic coordinates (>= 100000) in shift calculation
                                    if *svg_pos >= 100000 {
                                        continue;
                                    }
                                    let csv_pos = (*svg_pos as i32 + shift) as u32;
                                    if let Some((csv_base, _)) = pos_to_reactivity.get(&csv_pos) {
                                        total_checked += 1;
                                        // Compare normalized bases (T and U are equivalent)
                                        if normalize_base(*csv_base) == *base {
                                            match_count += 1;
                                        }
                                    }
                                }
                            }
                            
                            // Prefer shifts with more matches and more total checks
                            if total_checked > best_total_checked || 
                               (total_checked == best_total_checked && match_count > best_match_count) {
                                best_match_count = match_count;
                                best_total_checked = total_checked;
                                best_shift = shift;
                            }
                        }
                        
                        base_shifts.insert(*base, best_shift);
                        // DEBUG output removed to reduce noise
                        // eprintln!("Best shift for base {}: {} ({} matches out of {} checked, relative positions only)", 
                        //          base, best_shift, best_match_count, best_total_checked);
                    }
                }
    }

    // Create output file
    let output = File::create(output_file)?;
    let mut writer = BufWriter::new(output);

    // Extract template position to base mapping (most reliable source)
    // This uses "position.label in template: POS.BASE" format, which is the source of truth
    let template_pos2base = extract_template_pos_to_base(&svg_content)?;
    
    // Parse SVG template to extract circle positions with alignment
    let circles = parse_svg_template_for_circles(
        &svg_content,
        &pos_to_reactivity,
        &color_ranges,
        &config,
        &base_shifts,
        &template_pos2base,
    )?;
    // DEBUG output removed to reduce noise
    // println!("Found {} circles to draw", circles.len());
    // eprintln!(
    //     "DEBUG: circles_to_draw = {} (bases_filter = {})",
    //     circles.len(),
    //     bases_filter
    // );

    // Decide drawing order based on configuration
    // DEBUG output removed to reduce noise
    // eprintln!(
    //     "Drawing with circles_before_text: {}",
    //     config.circles_before_text
    // );
    if config.circles_before_text {
        // Draw color circles first, then skeleton and text
        // DEBUG output removed to reduce noise
        // eprintln!("Drawing circles first, then skeleton and text");
        draw_color_circles(&mut writer, &circles, &config)?;
        draw_svg_skeleton_and_text(&mut writer, &svg_content)?;
    } else {
        // Draw skeleton first (without text), then color circles, finally text
        // DEBUG output removed to reduce noise
        // eprintln!("Drawing skeleton first, then circles, then text");
        draw_svg_skeleton_without_text(&mut writer, &svg_content)?;
        draw_color_circles(&mut writer, &circles, &config)?;
        draw_svg_text_only(&mut writer, &svg_content, &config)?;
    }

    // Extract SVG width and height from template for legend positioning
    let svg_width = extract_svg_width(&svg_content);
    let svg_height = extract_svg_height(&svg_content);
    
    // Draw color bars (if enabled)
    if config.draw_color_bars {
        draw_color_bars(&mut writer, &base_groups, &color_ranges, bases_filter, &config, svg_width, svg_height)?;
    }

    // Write final closing SVG tag
    writeln!(writer, "</svg>")?;

    writer.flush()?;
    Ok(())
}

#[derive(Debug, Clone)]
struct SignalData {
    position: u32,
    stop_signal: f64,     // PIPC signal
    mutation_signal: f64, // RFC signal
    depth: usize,
}

#[derive(Debug, Clone)]
pub struct ReactivityData {
    pub position: u32,
    pub reactivity_stop: f64,
    pub reactivity_mutation: f64,
}

#[derive(Debug, Clone)]
struct GenePlotTask {
    gene_id: String,
    mod_signals: Vec<SignalData>,
    unmod_signals: Vec<SignalData>,
    reactivity_data: Option<Vec<ReactivityData>>,
    output_dir: String,
    is_high: bool,
}

#[derive(Debug, Clone)]
struct GeneStat {
    gene_id: String,
    coverage: f64,
    avg_depth: f64,
    cov_x_depth: f64,
}

fn load_signal_data(csv_file: &str) -> Result<HashMap<String, Vec<SignalData>>, Box<dyn Error>> {
    let mut gene_data: HashMap<String, Vec<SignalData>> = HashMap::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();

        if fields.len() >= 13 {
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].trim_matches('"').to_string();
            let position = fields[2].parse::<u32>()?;
            let rfc = fields[4].parse::<usize>()?;
            // fields[5] is pipe_truncation_count (float), needs to be converted to integer
            let pipc = fields[5].parse::<f64>()?.round() as usize;
            let depth = fields[6].parse::<usize>()?;

            // Calculate signal strength (normalized)
            let stop_signal = if depth > 0 {
                pipc as f64 / depth as f64
            } else {
                0.0
            };

            let mutation_signal = if depth > 0 {
                rfc as f64 / depth as f64
            } else {
                0.0
            };

            let signal_data = SignalData {
                position,
                stop_signal,
                mutation_signal,
                depth,
            };

            // Use chr_id+strand as gene_id to ensure separate processing of positive and negative strands
            let gene_id = format!("{}_{}", chr_id, strand);
            gene_data
                .entry(gene_id)
                .or_insert_with(Vec::new)
                .push(signal_data);
        }
    }

    // Sort data by position for each gene
    for data in gene_data.values_mut() {
        data.sort_by_key(|d| d.position);
    }

    Ok(gene_data)
}

pub fn load_reactivity_data(
    csv_file: &str,
) -> Result<HashMap<String, Vec<ReactivityData>>, Box<dyn Error>> {
    let mut gene_data: HashMap<String, Vec<ReactivityData>> = HashMap::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);

    // Determine file type (stop or mutation)
    let is_stop_file = csv_file.contains("stop");
    let is_mutation_file = csv_file.contains("mutation");

    for (line_num, line) in reader.lines().enumerate() {
        if line_num == 0 {
            continue; // Skip header row
        }

        let line = line?;
        let fields: Vec<&str> = line.split(',').collect();

        if fields.len() >= 6 {
            // Contains reactivity columns
            let chr_id = fields[0].trim_matches('"').to_string();
            let strand = fields[1].trim_matches('"').to_string();
            let position = fields[2].parse::<u32>()?;

            // Simplified format: column 4 is reactivity_stop, column 5 is reactivity_mutation
            let stop_val = parse_f64_allow_nan(fields[4]).ok_or_else(|| format!("Invalid reactivity_stop value: {}", fields[4]))?;
            let mut_val = parse_f64_allow_nan(fields[5]).ok_or_else(|| format!("Invalid reactivity_mutation value: {}", fields[5]))?;
            let (reactivity_stop, reactivity_mutation) = (stop_val, mut_val);

            let reactivity_data = ReactivityData {
                position,
                reactivity_stop,
                reactivity_mutation,
            };

            // Use chr_id+strand as gene_id to ensure separate processing of positive and negative strands
            let gene_id = format!("{}_{}", chr_id, strand);
            gene_data
                .entry(gene_id)
                .or_insert_with(Vec::new)
                .push(reactivity_data);
        }
    }

    // Sort data by position for each gene
    for data in gene_data.values_mut() {
        data.sort_by_key(|d| d.position);
    }

    Ok(gene_data)
}

fn gene_stats(gene_data: &HashMap<String, Vec<SignalData>>) -> HashMap<String, GeneStat> {
    let mut stats = HashMap::new();
    for (gene_id, signals) in gene_data {
        let n = signals.len();
        if n == 0 {
            continue;
        }
        let covered = signals.iter().filter(|d| d.depth > 0).count();
        let avg_depth = signals.iter().map(|d| d.depth as f64).sum::<f64>() / n as f64;
        let coverage = covered as f64 / n as f64;
        let cov_x_depth = coverage * avg_depth;
        stats.insert(
            gene_id.clone(),
            GeneStat {
                gene_id: gene_id.clone(),
                coverage,
                avg_depth,
                cov_x_depth,
            },
        );
    }
    stats
}

/// Parse GFF/GTF file, return gene_id->strand mapping
/// Supports GFF3 format: locus_tag=XXX or ID=XXX
/// Supports GTF format: gene_id "XXX" or transcript_id "XXX"
pub fn parse_gff_strand_map(gff_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut map = HashMap::new();
    let file = File::open(gff_file)?;
    let reader = BufReader::new(file);

    // Detect file format
    let mut is_gtf = false;
    let mut first_data_line = true;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Detect format (based on first data line)
        if first_data_line {
            let attr = fields[8];
            // GTF format: gene_id "ENSG00000243485"
            // GFF3 format: ID=BOGIMF_00001
            is_gtf = attr.contains("gene_id \"") || attr.contains("transcript_id \"");
            first_data_line = false;
        }

        let strand = fields[6].trim();
        let attr = fields[8];

        // Extract gene_id based on format
        let mut gene_id = None;

        if is_gtf {
            // GTF format parsing - prioritize transcript_id because norm files use transcript ID
            for kv in attr.split(';') {
                let kv = kv.trim();
                if kv.starts_with("transcript_id \"") {
                    // Prioritize transcript_id
                    if let Some(start) = kv.find('"') {
                        if let Some(end) = kv[start + 1..].find('"') {
                            gene_id = Some(kv[start + 1..start + 1 + end].to_string());
                            break;
                        }
                    }
                } else if kv.starts_with("gene_id \"") {
                    // If no transcript_id, use gene_id
                    if gene_id.is_none() {
                        if let Some(start) = kv.find('"') {
                            if let Some(end) = kv[start + 1..].find('"') {
                                gene_id = Some(kv[start + 1..start + 1 + end].to_string());
                            }
                        }
                    }
                }
            }
        } else {
            // GFF3 format parsing
            for kv in attr.split(';') {
                let kv = kv.trim();
                if kv.starts_with("locus_tag=") {
                    gene_id = Some(kv.trim_start_matches("locus_tag=").to_string());
                    break;
                } else if kv.starts_with("ID=") {
                    gene_id = Some(kv.trim_start_matches("ID=").to_string());
                }
            }
        }

        if let Some(id) = gene_id {
            map.insert(id, strand.to_string());
        }
    }

    Ok(map)
}

// Modified plot_signal_distributions function signature, added gff_file parameter
pub fn plot_signal_distributions(
    mod_file: &str,
    unmod_file: &str,
    output_dir: &str,
    num_threads: Option<usize>,
    coverage_threshold: Option<f64>,
    depth_threshold: Option<f64>,
    reactivity_file: Option<&str>,
    gff_file: Option<&str>, // New parameter
    logger: &mut crate::Logger,
) -> Result<(), Box<dyn Error>> {
    let start_time = std::time::Instant::now();

    // Record environment information and parameters
    logger.log("=== ModDetector Plot Function Log ===")?;
    logger.log(&format!("Software Version: v{}", crate::VERSION))?;
    logger.log(&format!(
        "Runtime: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S")
    ))?;
    logger.log(&format!("Modified Sample File: {}", mod_file))?;
    logger.log(&format!("Unmodified Sample File: {}", unmod_file))?;
    logger.log(&format!("Output Directory: {}", output_dir))?;
    if let Some(coverage) = coverage_threshold {
        logger.log(&format!("Coverage Threshold: {}", coverage))?;
    } else {
        logger.log("Coverage Threshold: Default 0.2")?;
    }
    if let Some(depth) = depth_threshold {
        logger.log(&format!("Depth Threshold: {}", depth))?;
    } else {
        logger.log("Depth Threshold: Default 50")?;
    }
    if let Some(reactivity_file) = reactivity_file {
        logger.log(&format!("Reactivity File: {}", reactivity_file))?;
    }
    if let Some(gff) = gff_file {
        logger.log(&format!("GFF Annotation File: {}", gff))?;
    }
    if let Some(threads) = num_threads {
        logger.log(&format!("Threads: {}", threads))?;
    } else {
        logger.log("Threads: Auto-detect")?;
    }
    logger.log("Starting signal distribution plotting...")?;

    // Display data loading information
    println!("[Loading data]");
    println!("    mod: {}", mod_file);
    println!("    unmod: {}", unmod_file);
    if let Some(reactivity_file) = reactivity_file {
        println!("    reactivity: {}", reactivity_file);
    }
    if let Some(gff_file) = gff_file {
        println!("    gene annotation: {}", gff_file);
    }
    println!();

    // Load data
    let mut mod_data = load_signal_data(mod_file)?;
    let mut unmod_data = load_signal_data(unmod_file)?;

    // Save unfiltered gene count
    let total_genes_before_filter = mod_data.len();
    let mut plus_strand_before_filter = 0;
    let mut minus_strand_before_filter = 0;
    for gene_id in mod_data.keys() {
        if gene_id.ends_with('+') {
            plus_strand_before_filter += 1;
        } else if gene_id.ends_with('-') {
            minus_strand_before_filter += 1;
        }
    }

    // ====== New: GFF annotation strand filtering ======
    if let Some(gff_path) = gff_file {
        let gff_map = parse_gff_strand_map(gff_path)?;
        // Only keep annotated strands
        let filter_keys: Vec<String> = mod_data
            .keys()
            .filter(|k| {
                // k format: geneid[_EupX_EdownY]_strand
                if let Some(pos) = k.rfind('_') {
                    let mut geneid = &k[..pos];
                    let strand = &k[pos + 1..];
                    // Remove _Eup/Edown and other upstream/downstream information
                    if let Some(eup_pos) = geneid.find("_Eup") {
                        geneid = &geneid[..eup_pos];
                    }
                    if let Some(anno_strand) = gff_map.get(geneid) {
                        return strand == anno_strand;
                    }
                }
                false
            })
            .cloned()
            .collect();
        mod_data = mod_data
            .into_iter()
            .filter(|(k, _)| filter_keys.contains(k))
            .collect();
        unmod_data = unmod_data
            .into_iter()
            .filter(|(k, _)| filter_keys.contains(k))
            .collect();
    }

    // Load reactivity data
    let mut reactivity_combined_data: Option<HashMap<String, Vec<ReactivityData>>> = None;

    if let Some(reactivity_file_path) = reactivity_file {
        if std::path::Path::new(reactivity_file_path).exists() {
            let reactivity_data = load_reactivity_data(reactivity_file_path)?;
            reactivity_combined_data = Some(reactivity_data);
        }
    }

    // Ensure output directory exists
    std::fs::create_dir_all(output_dir)?;

    // Set thread count
    let threads = num_threads.unwrap_or(rayon::current_num_threads());
    if let Some(threads) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap_or_else(|_| {
                println!("Warning: Unable to set thread count, using default settings")
            });
    }

    // Calculate gene coverage and depth
    let stats = gene_stats(&mod_data);

    // Set coverage and depth thresholds
    let coverage_threshold = coverage_threshold.unwrap_or(0.2); // Default 0.2
    let depth_threshold = depth_threshold.unwrap_or(50.0); // Default 50 reads

    // Display parameter information
    println!("[Params]");
    println!("    Threads: {}.", threads);
    println!(
        "    Coverage threshold: {:.1} %.",
        coverage_threshold * 100.0
    );
    println!("    Depth threshold: {:.0} reads.", depth_threshold);
    println!();

    // Count high/low group numbers
    let mut high_count = 0;
    let mut low_count = 0;
    for stat in stats.values() {
        if stat.coverage >= coverage_threshold && stat.avg_depth >= depth_threshold {
            high_count += 1;
        } else {
            low_count += 1;
        }
    }

    // Count positive/negative strand numbers (after filtering)
    let mut plus_strand = 0;
    let mut minus_strand = 0;
    for gene_id in mod_data.keys() {
        if gene_id.ends_with('+') {
            plus_strand += 1;
        } else if gene_id.ends_with('-') {
            minus_strand += 1;
        }
    }

    // Display data information
    println!("[Data info]");
    println!(
        "    Total gene regions: {}, {} (+ strand), {} (- strand).",
        total_genes_before_filter, plus_strand_before_filter, minus_strand_before_filter
    );
    if gff_file.is_some() {
        println!(
            "    Gene annotation matched: {}, {} (+ strand), {} (- strand).",
            mod_data.len(),
            plus_strand,
            minus_strand
        );
    }
    println!(
        "    High coverage: {}. Low coverage: {}.",
        high_count, low_count
    );
    println!();

    // Display progress information
    print!(
        "\r[Progressing] Plotting {}/{} ({:.1}%)",
        0,
        mod_data.len(),
        0.0
    );
    std::io::stdout().flush()?;

    // Create high/low directories
    let high_dir = format!("{}/high", output_dir);
    let low_dir = format!("{}/low", output_dir);
    std::fs::create_dir_all(&high_dir)?;
    std::fs::create_dir_all(&low_dir)?;

    // Prepare plotting tasks
    let mut plot_tasks: Vec<GenePlotTask> = Vec::new();
    for (gene_id, mod_signals) in &mod_data {
        if let Some(unmod_signals) = unmod_data.get(gene_id) {
            let stat = stats.get(gene_id).unwrap();
            let is_high = stat.coverage >= coverage_threshold && stat.avg_depth >= depth_threshold;

            // Get reactivity data
            let reactivity_data = if let Some(ref combined_data) = reactivity_combined_data {
                if let Some(reactivity) = combined_data.get(gene_id) {
                    Some(reactivity.clone())
                } else {
                    None
                }
            } else {
                None
            };

            plot_tasks.push(GenePlotTask {
                gene_id: gene_id.clone(),
                mod_signals: mod_signals.clone(),
                unmod_signals: unmod_signals.clone(),
                reactivity_data,
                output_dir: if is_high {
                    high_dir.clone()
                } else {
                    low_dir.clone()
                },
                is_high,
            });
        }
    }

    // Use atomic counter to track progress
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    let completed = Arc::new(AtomicUsize::new(0));
    let total_tasks = plot_tasks.len();

    // Process plotting tasks in parallel
    let results: Vec<Result<(), Box<dyn Error + Send + Sync>>> = plot_tasks
        .into_par_iter()
        .map(|task| {
            let gene_id = &task.gene_id;

            // Plot stop signal graph
            plot_stop_signals(
                gene_id,
                &task.mod_signals,
                &task.unmod_signals,
                &task.output_dir,
            )?;

            // Plot mutation signal graph
            plot_mutation_signals(
                gene_id,
                &task.mod_signals,
                &task.unmod_signals,
                &task.output_dir,
            )?;

            // Plot reactivity graph
            if let Some(ref reactivity_data) = task.reactivity_data {
                plot_reactivity_stop_signals(gene_id, reactivity_data, &task.output_dir)?;
                plot_reactivity_mutation_signals(gene_id, reactivity_data, &task.output_dir)?;
            }

            // Update progress
            let completed_count = completed.fetch_add(1, Ordering::Relaxed) + 1;
            let percentage = (completed_count as f64 * 100.0) / total_tasks as f64;
            print!(
                "\r[Progressing] Plotting {}/{} ({:.1}%)",
                completed_count, total_tasks, percentage
            );
            std::io::stdout().flush()?;

            Ok(())
        })
        .collect();

    // Check results
    let mut success_count = 0;
    let mut error_count = 0;
    for result in results {
        match result {
            Ok(_) => success_count += 1,
            Err(e) => {
                error_count += 1;
                eprintln!("Plotting error: {}", e);
            }
        }
    }

    let elapsed = start_time.elapsed();

    // Overwrite progress display, show output information
    println!(
        "\r[Outputs]  Success: {}. Failed: {}.",
        success_count, error_count
    );
    println!("    summary plot: {}/overall_scatter.png", output_dir);
    println!("    each gene plots: {}/", output_dir);
    println!("{}", crate::progress::format_time_used(elapsed));

    // Plot overall scatter
    let scatter_path = format!("{}/overall_scatter.png", output_dir);
    plot_overall_scatter(&stats, coverage_threshold, depth_threshold, &scatter_path)?;

    // Log completion
    logger.log(&format!("Signal distribution plotting completed"))?;
    logger.log(&format!("Successfully plotted genes: {}", success_count))?;
    if error_count > 0 {
        logger.log(&format!("Failed to plot genes: {}", error_count))?;
    }
    logger.log(&format!("Output Directory: {}", output_dir))?;
    logger.log(&format!("Total time: {:.2}s", elapsed.as_secs_f64()))?;

    if error_count > 0 {
        return Err(format!("{} gene regions failed to plot", error_count).into());
    }

    Ok(())
}

fn plot_stop_signals(
    gene_id: &str,
    mod_signals: &[SignalData],
    unmod_signals: &[SignalData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!(
        "{}/{}_stop_signals.png",
        output_dir,
        gene_id.replace('/', "_")
    );
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(10, 10, 10, 10);

    // Determine data range
    let all_positions: Vec<u32> = mod_signals
        .iter()
        .chain(unmod_signals.iter())
        .map(|d| d.position)
        .collect();

    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);

    let max_stop_signal = mod_signals
        .iter()
        .chain(unmod_signals.iter())
        .map(|d| d.stop_signal)
        .fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("{} - Stop Signal Distribution", gene_id),
            ("sans-serif", 30),
        )
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(*min_pos..*max_pos, 0.0..max_stop_signal * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Stop Signal (PIPC/Depth)")
        .draw()?;

    // Plot modified sample data (green)
    if !mod_signals.is_empty() {
        let mod_points: Vec<(u32, f64)> = mod_signals
            .iter()
            .map(|d| (d.position, d.stop_signal))
            .collect();

        chart
            .draw_series(LineSeries::new(
                mod_points.iter().map(|&(x, y)| (x, y)),
                GREEN.mix(0.8).stroke_width(2),
            ))?
            .label("Modified")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

        chart.draw_series(
            mod_points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 2, GREEN.mix(0.8).filled())),
        )?;
    }

    // Plot unmodified sample data (blue)
    if !unmod_signals.is_empty() {
        let unmod_points: Vec<(u32, f64)> = unmod_signals
            .iter()
            .map(|d| (d.position, d.stop_signal))
            .collect();

        chart
            .draw_series(LineSeries::new(
                unmod_points.iter().map(|&(x, y)| (x, y)),
                BLUE.mix(0.8).stroke_width(2),
            ))?
            .label("Unmodified")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

        chart.draw_series(
            unmod_points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 2, BLUE.mix(0.8).filled())),
        )?;
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

fn plot_mutation_signals(
    gene_id: &str,
    mod_signals: &[SignalData],
    unmod_signals: &[SignalData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!(
        "{}/{}_mutation_signals.png",
        output_dir,
        gene_id.replace('/', "_")
    );
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(10, 10, 10, 10);

    // Determine data range
    let all_positions: Vec<u32> = mod_signals
        .iter()
        .chain(unmod_signals.iter())
        .map(|d| d.position)
        .collect();

    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);

    let max_mutation_signal = mod_signals
        .iter()
        .chain(unmod_signals.iter())
        .map(|d| d.mutation_signal)
        .fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("{} - Mutation Signal Distribution", gene_id),
            ("sans-serif", 30),
        )
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(*min_pos..*max_pos, 0.0..max_mutation_signal * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Mutation Signal (RFC/Depth)")
        .draw()?;

    // Plot modified sample data (green)
    if !mod_signals.is_empty() {
        let mod_points: Vec<(u32, f64)> = mod_signals
            .iter()
            .map(|d| (d.position, d.mutation_signal))
            .collect();

        chart
            .draw_series(LineSeries::new(
                mod_points.iter().map(|&(x, y)| (x, y)),
                GREEN.mix(0.8).stroke_width(2),
            ))?
            .label("Modified")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

        chart.draw_series(
            mod_points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 2, GREEN.mix(0.8).filled())),
        )?;
    }

    // Plot unmodified sample data (blue)
    if !unmod_signals.is_empty() {
        let unmod_points: Vec<(u32, f64)> = unmod_signals
            .iter()
            .map(|d| (d.position, d.mutation_signal))
            .collect();

        chart
            .draw_series(LineSeries::new(
                unmod_points.iter().map(|&(x, y)| (x, y)),
                BLUE.mix(0.8).stroke_width(2),
            ))?
            .label("Unmodified")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

        chart.draw_series(
            unmod_points
                .iter()
                .map(|&(x, y)| Circle::new((x, y), 2, BLUE.mix(0.8).filled())),
        )?;
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

fn plot_reactivity_stop_signals(
    gene_id: &str,
    reactivity_data: &[ReactivityData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!(
        "{}/{}_reactivity_stop_signals.png",
        output_dir,
        gene_id.replace('/', "_")
    );
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(10, 10, 10, 10);

    // Determine data range
    let all_positions: Vec<u32> = reactivity_data.iter().map(|d| d.position).collect();
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);

    let min_pos_f = *min_pos as f64;
    let max_pos_f = *max_pos as f64;

    let max_reactivity = reactivity_data
        .iter()
        .map(|d| d.reactivity_stop)
        .fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("{} - Stop Reactivity Signal Distribution", gene_id),
            ("sans-serif", 30),
        )
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_pos_f..max_pos_f, 0.0..max_reactivity * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Stop Reactivity Signal")
        .draw()?;

    // Plot stop reactivity bar chart
    for data in reactivity_data {
        let x = data.position as f64;
        let y_stop = data.reactivity_stop;

        // Select color based on reactivity value
        let stop_color = if y_stop <= 0.2 {
            RGBColor(128, 128, 128) // Gray
        } else if y_stop <= 0.4 {
            RGBColor(255, 255, 0) // Yellow
        } else {
            RGBColor(255, 0, 0) // Red
        };

        // Plot stop reactivity bar chart
        if y_stop > 0.0 {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x - 0.2, 0.0), (x + 0.2, y_stop)],
                stop_color.filled(),
            )))?;
        }
    }

    // Add legend
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(128, 128, 128).filled(),
        )))?
        .label("0-0.2 (Gray)");

    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(255, 255, 0).filled(),
        )))?
        .label("0.2-0.4 (Yellow)");

    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(255, 0, 0).filled(),
        )))?
        .label("0.4-1.0 (Red)");

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

fn plot_reactivity_mutation_signals(
    gene_id: &str,
    reactivity_data: &[ReactivityData],
    output_dir: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let filename = format!(
        "{}/{}_reactivity_mutation_signals.png",
        output_dir,
        gene_id.replace('/', "_")
    );
    let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let root = root.margin(10, 10, 10, 10);

    // Determine data range
    let all_positions: Vec<u32> = reactivity_data.iter().map(|d| d.position).collect();
    let min_pos = all_positions.iter().min().unwrap_or(&0);
    let max_pos = all_positions.iter().max().unwrap_or(&0);

    let min_pos_f = *min_pos as f64;
    let max_pos_f = *max_pos as f64;

    let max_reactivity = reactivity_data
        .iter()
        .map(|d| d.reactivity_mutation)
        .fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("{} - Mutation Reactivity Signal Distribution", gene_id),
            ("sans-serif", 30),
        )
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_pos_f..max_pos_f, 0.0..max_reactivity * 1.1)?;

    chart
        .configure_mesh()
        .x_desc("Position")
        .y_desc("Mutation Reactivity Signal")
        .draw()?;

    // Plot mutation reactivity bar chart
    for data in reactivity_data {
        let x = data.position as f64;
        let y_mut = data.reactivity_mutation;

        // Select color based on reactivity value
        let mut_color = if y_mut <= 0.2 {
            RGBColor(128, 128, 128) // Gray
        } else if y_mut <= 0.4 {
            RGBColor(255, 255, 0) // Yellow
        } else {
            RGBColor(255, 0, 0) // Red
        };

        // Plot mutation reactivity bar chart
        if y_mut > 0.0 {
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x - 0.2, 0.0), (x + 0.2, y_mut)],
                mut_color.filled(),
            )))?;
        }
    }

    // Add legend
    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(128, 128, 128).filled(),
        )))?
        .label("0-0.2 (Gray)");

    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(255, 255, 0).filled(),
        )))?
        .label("0.2-0.4 (Yellow)");

    chart
        .draw_series(std::iter::once(Rectangle::new(
            [(0.0, 0.0), (10.0, 0.0)],
            RGBColor(255, 0, 0).filled(),
        )))?
        .label("0.4-1.0 (Red)");

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}

fn plot_overall_scatter(
    stats: &HashMap<String, GeneStat>,
    coverage_threshold: f64,
    depth_threshold: f64,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(output_path, (1200, 900)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut all_coverage: Vec<f64> = Vec::new();
    let mut all_depth: Vec<f64> = Vec::new();
    for s in stats.values() {
        all_coverage.push(s.coverage);
        all_depth.push(s.avg_depth);
    }

    // coverage does not use log10 transformation, depth still uses log10 transformation
    let log_depth: Vec<f64> = all_depth.iter().map(|&x| (x + 1.0).log10()).collect();

    // Force x-axis and y-axis to start from 0
    let min_cov = 0.0;
    let max_cov = *all_coverage
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(&1.0);
    let min_depth = 0.0;
    let max_depth = *log_depth
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(&2.0);

    let mut chart = ChartBuilder::on(&root)
        .caption("Gene Coverage-Depth Scatter", ("sans-serif", 36))
        .x_label_area_size(60)
        .y_label_area_size(80)
        .margin(20)
        .build_cartesian_2d(min_cov..max_cov, min_depth..max_depth)?;

    chart
        .configure_mesh()
        .x_desc("Coverage")
        .y_desc("log10(Average Depth)")
        .draw()?;

    // Draw threshold lines (coverage does not use log10, depth uses log10)
    let log_depth_threshold = (depth_threshold + 1.0).log10();

    // Vertical coverage threshold line
    chart.draw_series(LineSeries::new(
        vec![
            (coverage_threshold, min_depth),
            (coverage_threshold, max_depth),
        ],
        RGBColor(255, 0, 0).mix(0.7).stroke_width(2),
    ))?;

    // Horizontal depth threshold line
    chart.draw_series(LineSeries::new(
        vec![
            (min_cov, log_depth_threshold),
            (max_cov, log_depth_threshold),
        ],
        RGBColor(0, 0, 255).mix(0.7).stroke_width(2),
    ))?;

    // Count high/low group numbers
    let mut high_count = 0;
    let mut low_count = 0;
    for s in stats.values() {
        if s.coverage >= coverage_threshold && s.avg_depth >= depth_threshold {
            high_count += 1;
        } else {
            low_count += 1;
        }
    }

    // Draw points
    for s in stats.values() {
        let cov = s.coverage;
        let log_depth = (s.avg_depth + 1.0).log10();
        let color = if s.coverage >= coverage_threshold && s.avg_depth >= depth_threshold {
            RED
        } else {
            RGBColor(180, 180, 180)
        };
        chart.draw_series(std::iter::once(Circle::new(
            (cov, log_depth),
            7,
            color.filled(),
        )))?;
    }

    // Add legend and count annotations
    chart
        .draw_series(std::iter::once(Circle::new((0.0, 0.0), 0, RED.filled())))?
        .label(&format!("High (N={})", high_count));
    chart
        .draw_series(std::iter::once(Circle::new(
            (0.0, 0.0),
            0,
            RGBColor(180, 180, 180).filled(),
        )))?
        .label(&format!("Low (N={})", low_count));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    // Add threshold annotations
    chart.draw_series(std::iter::once(Text::new(
        format!("Coverage ≥ {:.3}", coverage_threshold),
        (max_cov - 0.1, max_depth - 0.3),
        ("sans-serif", 20).into_font().color(&RGBColor(255, 0, 0)),
    )))?;

    chart.draw_series(std::iter::once(Text::new(
        format!("Depth ≥ {:.3}", depth_threshold),
        (min_cov + 0.1, max_depth - 0.3),
        ("sans-serif", 20).into_font().color(&RGBColor(0, 0, 255)),
    )))?;

    root.present()?;
    Ok(())
}

// ==================== SVG Plotting Functions ====================

#[derive(Debug, Clone)]
pub struct SvgReactivityData {
    pub position: u32,
    pub base: char,
    pub reactivity: f64,
}

/// Base metadata for metadata-driven visualization architecture
/// This structure contains all information needed for precise, atomic updates
#[derive(Debug, Clone)]
pub struct BaseMetadata {
    pub id: u32,               // Position ID (unique identifier)
    pub base: char,            // Normalized base type (A, C, G, U - T is normalized to U)
    pub x: f64,                // SVG x coordinate
    pub y: f64,                // SVG y coordinate
    pub reactivity: f64,       // Original reactivity value
}

#[derive(Debug, Clone)]
pub struct ColorRange {
    pub min: f64,
    pub max: f64,
    pub color: String,
}

#[derive(Debug, Clone)]
pub struct BaseColorRanges {
    pub a: Vec<ColorRange>,
    pub t: Vec<ColorRange>,
    pub c: Vec<ColorRange>,
    pub g: Vec<ColorRange>,
}

impl Default for BaseColorRanges {
    fn default() -> Self {
        Self {
            a: vec![
                ColorRange {
                    min: 0.0,
                    max: 0.3,
                    color: "rgba(245, 245, 245, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.3,
                    max: 0.65,
                    color: "rgba(255, 165, 0, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.65,
                    max: 1.0,
                    color: "rgba(255, 0, 0, 0.8)".to_string(),
                },
            ],
            t: vec![
                ColorRange {
                    min: 0.0,
                    max: 0.2,
                    color: "rgba(245, 245, 245, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.2,
                    max: 0.55,
                    color: "rgba(255, 165, 0, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.55,
                    max: 1.0,
                    color: "rgba(255, 0, 0, 0.8)".to_string(),
                },
            ],
            c: vec![
                ColorRange {
                    min: 0.0,
                    max: 0.3,
                    color: "rgba(245, 245, 245, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.3,
                    max: 0.6,
                    color: "rgba(255, 165, 0, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.6,
                    max: 1.0,
                    color: "rgba(255, 0, 0, 0.8)".to_string(),
                },
            ],
            g: vec![
                ColorRange {
                    min: 0.0,
                    max: 0.2,
                    color: "rgba(245, 245, 245, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.2,
                    max: 0.5,
                    color: "rgba(255, 165, 0, 0.8)".to_string(),
                },
                ColorRange {
                    min: 0.5,
                    max: 1.0,
                    color: "rgba(255, 0, 0, 0.8)".to_string(),
                },
            ],
        }
    }
}

/// Load reactivity data from CSV file for SVG plotting
pub fn load_svg_reactivity_data(
    csv_file: &str,
    signal_type: &str,
    strand_filter: &str,
) -> Result<Vec<SvgReactivityData>, Box<dyn Error>> {
    let mut data = Vec::new();
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<Vec<_>, _>>()?;

    if lines.is_empty() {
        return Ok(data);
    }

    // Auto-detect CSV format from first line
    let first_line = &lines[0];
    let fields: Vec<&str> = first_line.split(',').collect();

    let has_header = if fields.len() >= 3 {
        fields[2].parse::<u32>().is_err() // If position field is not numeric, it's a header
    } else {
        false
    };

    let has_strand_column = fields.len() >= 5;

    // Auto-detect signal columns based on header
    let mut signal_columns: Vec<(String, usize)> = Vec::new();
    if has_header {
        // Parse header to find signal columns
        for (i, field) in fields.iter().enumerate() {
            if field.contains("reactivity_") {
                let signal_name = field.replace("reactivity_", "");
                signal_columns.push((signal_name, i));
            } else if i >= 4 && *field != "Strand" && *field != "Position" && *field != "Base" {
                // Assume this is a signal column if it's not a standard column
                signal_columns.push((field.to_string(), i));
            }
        }
    } else {
        // For files without header, assume last column is the signal
        if fields.len() >= 5 {
            signal_columns.push(("score".to_string(), fields.len() - 1));
        }
    }

    // If no signal columns detected, use default
    if signal_columns.is_empty() {
        signal_columns.push((
            "score".to_string(),
            if fields.len() == 6 {
                4
            } else if fields.len() == 5 {
                4
            } else {
                3
            },
        ));
    }

    // Check if user-requested column name exists; if not and only one signal column exists, allow fallback to that unique column
    let requested_exists =
        signal_type == "all" || signal_columns.iter().any(|(name, _)| name == signal_type);

    // Process data lines
    let start_line = if has_header { 1 } else { 0 };

    for line in &lines[start_line..] {
        let fields: Vec<&str> = line.split(',').collect();

        // Skip empty lines
        if fields.is_empty() || fields[0].trim().is_empty() {
            continue;
        }

        // Check strand filter
        let should_process = if has_strand_column && fields.len() >= 5 {
            match strand_filter {
                "+" => fields[1].trim() == "+",
                "-" => fields[1].trim() == "-",
                "both" => true,
                _ => fields[1].trim() == "+", // Default to positive strand
            }
        } else {
            true // No strand column, process all data
        };

        if should_process {
            // Process each signal type
            for (signal_name, column_index) in &signal_columns {
                // Check if this signal type should be processed
                let should_process_signal = if signal_type == "all" {
                    true
                } else if signal_name == signal_type {
                    true
                } else if !requested_exists && signal_columns.len() == 1 {
                    // Auto fallback: when user-specified column is not found and only one signal column exists, use that unique column
                    true
                } else {
                    false
                };

                if should_process_signal && column_index < &fields.len() {
                    if let Ok(position) = fields[2].parse::<u32>() {
                        if let Some(reactivity) = parse_f64_allow_nan(fields[*column_index]) {
                            let base = normalize_base(fields[3].chars().next().unwrap_or('N'));
                            data.push(SvgReactivityData {
                                position,
                                base,
                                reactivity,
                            });
                        }
                    }
                }
            }
        }
    }

    Ok(data)
}

/// Normalize base character (U -> T, uppercase)
fn normalize_base(base: char) -> char {
    let base = base.to_uppercase().next().unwrap_or('N');
    if base == 'U' {
        'T'
    } else {
        base
    }
}

/// Extract template position to base mapping from SVG template
/// Uses the most reliable format: "position.label in template: POS.BASE"
/// This is the source of truth for base types in the SVG template
fn extract_template_pos_to_base(svg_content: &str) -> Result<HashMap<u32, char>, Box<dyn Error>> {
    let mut pos_to_base = HashMap::new();
    // Match format: "position.label in template: POS.BASE"
    // Example: "position.label in template: 1792.G"
    let re = Regex::new(r"position\.label in template:\s*(\d+)\.([ATCGU])")?;
    
    for cap in re.captures_iter(svg_content) {
        let pos: u32 = cap[1].parse()?;
        let base_char = cap[2].chars().next().ok_or("Invalid base character")?;
        let normalized_base = normalize_base(base_char);
        pos_to_base.insert(pos, normalized_base);
    }
    
    Ok(pos_to_base)
}

/// Extract SVG width from template (from viewBox or width attribute)
fn extract_svg_width(svg_content: &str) -> Option<f64> {
    // Try to extract from viewBox first (format: viewBox="x y width height")
    let viewbox_regex = Regex::new(r#"viewBox\s*=\s*"[\s]*[\d.]+[\s]+[\d.]+[\s]+([\d.]+)[\s]+[\d.]+""#).ok()?;
    if let Some(caps) = viewbox_regex.captures(svg_content) {
        if let Ok(width) = caps[1].parse::<f64>() {
            return Some(width);
        }
    }
    
    // Try to extract from width attribute (format: width="xxx" or width='xxx')
    let width_regex = Regex::new(r#"width\s*=\s*["']([\d.]+)["']"#).ok()?;
    if let Some(caps) = width_regex.captures(svg_content) {
        if let Ok(width) = caps[1].parse::<f64>() {
            return Some(width);
        }
    }
    
    None
}

/// Extract SVG height from template (from viewBox or height attribute)
fn extract_svg_height(svg_content: &str) -> Option<f64> {
    // Try to extract from viewBox first (format: viewBox="x y width height")
    let viewbox_regex = Regex::new(r#"viewBox\s*=\s*"[\s]*[\d.]+[\s]+[\d.]+[\s]+[\d.]+[\s]+([\d.]+)""#).ok()?;
    if let Some(caps) = viewbox_regex.captures(svg_content) {
        if let Ok(height) = caps[1].parse::<f64>() {
            return Some(height);
        }
    }
    
    // Try to extract from height attribute (format: height="xxx" or height='xxx')
    let height_regex = Regex::new(r#"height\s*=\s*["']([\d.]+)["']"#).ok()?;
    if let Some(caps) = height_regex.captures(svg_content) {
        if let Ok(height) = caps[1].parse::<f64>() {
            return Some(height);
        }
    }
    
    None
}

/// Extract positions and bases from SVG template
fn extract_svg_positions_and_bases(svg_content: &str) -> Result<HashMap<u32, char>, Box<dyn Error>> {
    let position_regex = Regex::new(r"<title>(\d+)")?;
    // Priority: extract position from "position.label in template: XXX" if available
    let template_position_regex = Regex::new(r"position\.label in template: (\d+)")?;
    let base_regex = Regex::new(r"<text[^>]*>([ATCGU])</text>")?;
    let mut pos2base = HashMap::new();
    
    let lines: Vec<&str> = svg_content.lines().collect();
    
    for (idx, line) in lines.iter().enumerate() {
        if line.contains("title>") && !line.contains("numbering-label") {
            // Priority: use template position if available, otherwise use title position
            let svg_position_opt = if let Some(template_caps) = template_position_regex.captures(line) {
                template_caps[1].parse::<u32>().ok()
            } else if let Some(title_caps) = position_regex.captures(line) {
                title_caps[1].parse::<u32>().ok()
            } else {
                None
            };
            
            if let Some(svg_position) = svg_position_opt {
                // Try to find base in current line or nearby lines
                let mut base: Option<char> = None;
                
                // Check current line
                if let Some(base_caps) = base_regex.captures(line) {
                    base = base_caps[1].chars().next();
                } else {
                    // Check nearby lines (within 3 lines)
                    for i in (idx.saturating_sub(1))..(idx + 4).min(lines.len()) {
                        if let Some(base_caps) = base_regex.captures(lines[i]) {
                            base = base_caps[1].chars().next();
                            break;
                        }
                    }
                }
                
                if let Some(b) = base {
                    pos2base.insert(svg_position, normalize_base(b));
                }
            }
        }
    }
    
    Ok(pos2base)
}

/// Load reference sequence from FASTA or plain text file
pub fn load_reference_sequence(ref_path: &str) -> Result<HashMap<u32, char>, Box<dyn Error>> {
    let mut seq = String::new();
    let lower_path = ref_path.to_lowercase();

    if lower_path.ends_with(".fa")
        || lower_path.ends_with(".fasta")
        || lower_path.ends_with(".fa.gz")
        || lower_path.ends_with(".fasta.gz")
    {
        // Parse FASTA file
        let file = File::open(ref_path)?;
        let reader = BufReader::new(file);
        let mut in_sequence = false;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                in_sequence = true;
                continue;
            }
            if in_sequence && !line.trim().is_empty() {
                seq.push_str(&line.trim());
            }
        }
    } else {
        // Treat as plain text sequence
        let file = File::open(ref_path)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            if !line.trim().is_empty() && !line.starts_with('>') {
                seq.push_str(&line.trim());
            }
        }
    }

    if seq.is_empty() {
        return Err("Reference sequence is empty or not found".into());
    }

    // Convert to position->base map (1-based positions)
    let mut pos_to_base = HashMap::new();
    for (i, base) in seq.chars().enumerate() {
        pos_to_base.insert((i + 1) as u32, normalize_base(base));
    }

    Ok(pos_to_base)
}

/// Find best shift for a given base type by maximizing matches with reference sequence
pub fn find_best_shift_for_base(
    base_data: &[&SvgReactivityData],
    ref_pos2base: &HashMap<u32, char>,
    max_shift: u32,
) -> i32 {
    if base_data.is_empty() {
        return 0;
    }

    // Build observed map
    let obs: HashMap<u32, char> = base_data
        .iter()
        .map(|data| (data.position, data.base))
        .collect();

    let mut best_shift = 0;
    let mut best_match = -1;

    for s in -(max_shift as i32)..=(max_shift as i32) {
        let mut match_count = 0;

        // Compare positions that exist in both after shifting
        for (&pos, &base) in &obs {
            let ref_pos = (pos as i32 + s) as u32;
            if let Some(&ref_base) = ref_pos2base.get(&ref_pos) {
                if base == ref_base {
                    match_count += 1;
                }
            }
        }

        if match_count > best_match {
            best_match = match_count;
            best_shift = s;
        }
    }

    best_shift
}

/// Determine color based on reactivity value and color ranges
fn determine_color(shape_value: f64, color_ranges: &[ColorRange]) -> String {
    // Handle special case: -1.0 means no data
    if shape_value == -1.0 {
        return "rgba(169, 169, 169, 0.8)".to_string();
    }

    // Handle NaN or invalid values
    if shape_value.is_nan() || shape_value.is_infinite() {
        return "rgba(169, 169, 169, 0.8)".to_string();
    }

    // Check each range, using half-open intervals [min, max) except for the last range
    for (idx, range) in color_ranges.iter().enumerate() {
        let is_last_range = idx == color_ranges.len() - 1;
        
        if is_last_range {
            // Last range: include max boundary [min, max]
            if shape_value >= range.min && shape_value <= range.max {
                return range.color.clone();
            }
        } else {
            // Other ranges: half-open interval [min, max)
            if shape_value >= range.min && shape_value < range.max {
                return range.color.clone();
            }
        }
    }

    // If value is outside all ranges, map to the last range's color if > max, or first range if < min
    if !color_ranges.is_empty() {
        let first_range = &color_ranges[0];
        let last_range = &color_ranges[color_ranges.len() - 1];
        
        if shape_value < first_range.min {
            return first_range.color.clone();
        } else if shape_value > last_range.max {
            return last_range.color.clone();
        }
    }

    // Fallback: gray for no data
    "rgba(169, 169, 169, 0.8)".to_string()
}

/// Draw color bar legend for a base type
/// Uses structured SVG groups for easier manipulation in JavaScript
fn draw_color_bar(
    writer: &mut BufWriter<File>,
    color_ranges: &[ColorRange],
    base: char,
    position: usize,
    config: &SvgDrawConfig,
    svg_width: Option<f64>,
    svg_height: Option<f64>,
    total_bases: usize,
) -> Result<(), Box<dyn Error>> {
    // Use configurable dimensions for horizontal/vertical layouts
    let item_width = config.legend_item_width;  // Width of each color segment (for horizontal layout)
    let item_height = config.legend_item_height; // Height of each color segment (for vertical layout)
    let text_font_size = 8.0;
    let text_spacing = 3.0; // Spacing between text and rect
    let row_spacing = 20.0; // Spacing between rows
    let row_height = item_height + text_font_size + text_spacing + row_spacing;
    
    // Calculate legend width: -1 item + all color range items
    let legend_width = item_width * (1.0 + color_ranges.len() as f64);
    let margin = 20.0; // Margin from edge
    
    // Calculate initial position for horizontal layout (default)
    // Place legend on the right side to avoid blocking RNA structure
    // If SVG width is available, place legend on the right; otherwise use left side as fallback
    let base_x = if let Some(width) = svg_width {
        // Place on right side: SVG width - legend width - margin
        (width - legend_width - margin).max(margin)
    } else {
        // Fallback: use left side if width not available
        0.0
    };
    
    // Calculate y position from bottom
    // Place legend at bottom right to avoid blocking RNA structure
    let base_y = if let Some(height) = svg_height {
        // Calculate total legend height
        let total_legend_height = (total_bases as f64) * row_height;
        // Place from bottom: SVG height - total legend height - margin
        // Then add position offset for this specific base
        (height - total_legend_height - margin) + (position as f64) * row_height
    } else {
        // Fallback: use top if height not available
        (position as f64) * row_height + 50.0
    };
    
    // Create main legend group with class for easy selection and initial transform
    writeln!(writer, r#"<g id="legend-{}" class="legend-group" data-base="{}" transform="translate({}, {})">"#, base, base, base_x, base_y)?;
    
    // Draw base label (outside the items group for easier positioning)
    writeln!(
        writer,
        r#"<text x="0" y="{}" font-family="Arial" font-size="12" fill="black" font-weight="bold" class="legend-base-label">{}</text>"#,
        text_font_size, base
    )?;
    
    // Create items group - all color items will be in this group
    writeln!(writer, r#"<g class="legend-items" data-base="{}">"#, base)?;
    
    // Draw -1 value item (structured as a group)
    writeln!(writer, r#"<g class="legend-item" data-value="-1">"#)?;
    writeln!(
        writer,
        r#"<rect x="0" y="{}" width="{}" height="{}" fill="rgba(169, 169, 169, 0.8)" stroke="black" stroke-width="1" class="legend-rect" />"#,
        text_font_size + text_spacing, item_width, item_height
    )?;
    writeln!(
        writer,
        r#"<text x="{}" y="{}" font-family="Arial" font-size="{}" fill="black" class="legend-text" text-anchor="middle">-1</text>"#,
        item_width / 2.0,
        text_font_size,
        text_font_size
    )?;
    writeln!(writer, r#"</g>"#)?;
    
    // Draw color range items (each as a structured group)
    let mut current_x = item_width; // Start after -1 item
    for range in color_ranges.iter() {
        writeln!(writer, r#"<g class="legend-item" data-min="{:.2}" data-max="{:.2}">"#, range.min, range.max)?;
        // Rect positioned relative to item group
        writeln!(
            writer,
            r#"<rect x="{}" y="{}" width="{}" height="{}" fill="{}" stroke="black" stroke-width="1" class="legend-rect" />"#,
            current_x,
            text_font_size + text_spacing,
            item_width,
            item_height,
            range.color
        )?;
        // Text above rect, centered
        writeln!(
            writer,
            r#"<text x="{}" y="{}" font-family="Arial" font-size="{}" fill="black" class="legend-text" text-anchor="middle">{:.2} - {:.2}</text>"#,
            current_x + item_width / 2.0,
            text_font_size,
            text_font_size,
            range.min,
            range.max
        )?;
        writeln!(writer, r#"</g>"#)?;
        current_x += item_width; // Next item connects directly (no gap)
    }
    
    writeln!(writer, r#"</g>"#)?; // Close legend-items group
    writeln!(writer, r#"</g>"#)?; // Close legend group
    
    Ok(())
}

/// Plot reactivity data to SVG without alignment
pub fn plot_reactivity_to_svg(
    reactivity_file: &str,
    svg_template: &str,
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    println!("=== plot_reactivity_to_svg called ===");
    println!("reactivity_file: {}", reactivity_file);
    println!("signal_type: {}", signal_type);
    println!("strand_filter: {}", strand_filter);

    // Use default configuration, text before color circles (top layer)
    let config = SvgDrawConfig::default();

    // Call new drawing function (no reference sequence alignment)
    draw_svg_with_custom_layers(
        reactivity_file,
        svg_template,
        output_file,
        bases_filter,
        signal_type,
        strand_filter,
        color_ranges,
        config,
        None, // No reference sequence alignment
        5,    // Default max_shift
    )
}

/// Plot multiple signal types to SVG (auto-detect signals from CSV)
pub fn plot_multiple_signals_to_svg(
    reactivity_file: &str,
    svg_template: &str,
    output_dir: &str,
    bases_filter: &str,
    strand_filter: &str,
    ref_sequence_file: Option<&str>,
    max_shift: u32,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    // Auto-detect signal types from CSV file
    let signal_types = detect_signal_types(reactivity_file)?;

    // Create output directory
    std::fs::create_dir_all(output_dir)?;

    // Plot each signal type
    for signal_type in &signal_types {
        let output_file = format!("{}/rna_structure_colored_{}.svg", output_dir, signal_type);

        if let Some(ref_file) = ref_sequence_file {
            // Use reference sequence alignment
            plot_reactivity_to_svg_with_alignment(
                reactivity_file,
                svg_template,
                &output_file,
                bases_filter,
                signal_type,
                strand_filter,
                ref_file,
                max_shift,
                color_ranges.clone(),
            )?;
        } else {
            // No alignment
            plot_reactivity_to_svg(
                reactivity_file,
                svg_template,
                &output_file,
                bases_filter,
                signal_type,
                strand_filter,
                color_ranges.clone(),
            )?;
        }
    }

    Ok(())
}

/// Auto-detect signal types from CSV file header
fn detect_signal_types(csv_file: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let file = File::open(csv_file)?;
    let reader = BufReader::new(file);
    let first_line = reader.lines().next().ok_or("Empty CSV file")??;
    let fields: Vec<&str> = first_line.split(',').collect();

    let mut signal_types = Vec::new();
    for field in &fields {
        if field.contains("reactivity_") {
            let signal_name = field.replace("reactivity_", "");
            signal_types.push(signal_name);
        }
    }

    // If no reactivity columns found, assume single signal
    if signal_types.is_empty() {
        signal_types.push("score".to_string());
    }

    Ok(signal_types)
}

/// Plot reactivity data to SVG with alignment
pub fn plot_reactivity_to_svg_with_alignment(
    reactivity_file: &str,
    svg_template: &str,
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    ref_sequence_file: &str,
    max_shift: u32,
    color_ranges: Option<BaseColorRanges>,
) -> Result<(), Box<dyn Error>> {
    // Use default configuration, color circles before text
    let config = SvgDrawConfig::default();

    // Use reference sequence alignment
    draw_svg_with_custom_layers(
        reactivity_file,
        svg_template,
        output_file,
        bases_filter,
        signal_type,
        strand_filter,
        color_ranges,
        config,
        Some(ref_sequence_file),
        max_shift,
    )
}

// ============================================================================
// Interactive HTML Visualization Functions
// ============================================================================

/// Generate HTML for a single base control panel
fn generate_base_control_panel(base: char, name: &str, val1: f64, val2: f64) -> String {
    let val1_str = format!("{:.2}", val1);
    let val2_str = format!("{:.2}", val2);
    let base_str = base.to_string();
    format!(
        r##"
                <div class="base-control-panel" id="baseControl{}">
                    <div class="base-header">
                        <label><input type="checkbox" id="baseEnable{}" checked> {} ({})</label>
                    </div>
                    <div class="color-range-controls">
                        <div class="range-control">
                            <div class="distribution-label">Reactivity Distribution</div>
                            <div class="reactivity-distribution-preview" id="distributionPreview{}">
                                <canvas id="distributionCanvas{}"></canvas>
                            </div>
                            <label>Low-Mid Threshold</label>
                            <input type="range" id="range{}1" min="0" max="1" step="0.01" value="{}">
                            <span id="range{}1Value">{}</span>
                        </div>
                        <div class="range-control">
                            <label>Mid-High Threshold</label>
                            <input type="range" id="range{}2" min="0" max="1" step="0.01" value="{}">
                            <span id="range{}2Value">{}</span>
                        </div>
                        <div class="color-picker-group">
                            <div class="color-picker-item">
                                <label>Low Color:</label>
                                <input type="text" id="color{}1" value="#f5f5f5" class="spectrum-color-picker">
                                <span class="color-preview" id="color{}1Preview" style="background-color: #f5f5f5;"></span>
                            </div>
                            <div class="color-picker-item">
                                <label>Mid Color:</label>
                                <input type="text" id="color{}2" value="#ffa500" class="spectrum-color-picker">
                                <span class="color-preview" id="color{}2Preview" style="background-color: #ffa500;"></span>
                            </div>
                            <div class="color-picker-item">
                                <label>High Color:</label>
                                <input type="text" id="color{}3" value="#ff0000" class="spectrum-color-picker">
                                <span class="color-preview" id="color{}3Preview" style="background-color: #ff0000;"></span>
                            </div>
                        </div>
                    </div>
                </div>"##,
        base_str, base_str, name, base_str, base_str, base_str, base_str, val1_str, base_str, val1_str,
        base_str, val2_str, base_str, val2_str, base_str, base_str, base_str, base_str, base_str, base_str
    )
}

/// Generate an interactive HTML visualization with embedded SVG and JavaScript
pub fn generate_interactive_html(
    svg_content: &str,
    reactivity_data: &[SvgReactivityData],
    output_file: &str,
    signal_type: &str,
    color_ranges: &BaseColorRanges,
    initial_circle_filled: bool,
    base_metadata: Option<&[BaseMetadata]>,
) -> Result<(), Box<dyn Error>> {
    // Get version from crate::VERSION (defined in main.rs)
    let version = crate::VERSION;
    
    // Extract SVG filename from output file path for loading message
    let svg_filename = std::path::Path::new(output_file)
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| {
            // Try to extract meaningful name from filename
            // e.g., "br_k_m2018_mut_shapemap_EC_1M7_rep2_EC_16S.html" -> "EC_16S.svg"
            if let Some(gene_match) = s.rsplit('_').next() {
                if gene_match.starts_with("EC_") || gene_match.starts_with("Human_") {
                    format!("{}.svg", gene_match)
                } else {
                    format!("{}.svg", s)
                }
            } else {
                format!("{}.svg", s)
            }
        })
        .unwrap_or_else(|| "RNA structure.svg".to_string());
    
    let mut html = String::new();
    
    // HTML header with embedded CSS and JavaScript
    html.push_str(r###"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Interactive RNA Structure Visualization</title>
    <style>
        /* Initial page styles */
        #initialPage {
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
            z-index: 10000;
            transition: opacity 0.5s ease-out;
        }
        #initialPage.hidden {
            opacity: 0;
            pointer-events: none;
        }
        .initial-container {
            background: white;
            border-radius: 12px;
            padding: 40px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            max-width: 500px;
            width: 90%;
            text-align: center;
        }
        .initial-container h1 {
            color: #333;
            margin-bottom: 10px;
            font-size: 28px;
        }
        .initial-container p {
            color: #666;
            margin-bottom: 30px;
            font-size: 14px;
        }
        .file-upload-area {
            border: 2px dashed #667eea;
            border-radius: 8px;
            padding: 40px 20px;
            margin-bottom: 20px;
            cursor: pointer;
            transition: all 0.3s ease;
            background: #f8f9ff;
        }
        .file-upload-area:hover {
            border-color: #764ba2;
            background: #f0f2ff;
        }
        .file-upload-area.dragover {
            border-color: #764ba2;
            background: #e8ebff;
            transform: scale(1.02);
        }
        .file-upload-icon {
            font-size: 48px;
            color: #667eea;
            margin-bottom: 15px;
        }
        .file-upload-text {
            color: #333;
            font-size: 16px;
            margin-bottom: 10px;
        }
        .file-upload-hint {
            color: #999;
            font-size: 12px;
        }
        #fileInput {
            display: none;
        }
        .upload-button {
            background: #667eea;
            color: white;
            border: none;
            padding: 12px 30px;
            border-radius: 6px;
            font-size: 16px;
            cursor: pointer;
            transition: background 0.3s ease;
            margin-top: 10px;
        }
        .upload-button:hover {
            background: #764ba2;
        }
        .upload-button:disabled {
            background: #ccc;
            cursor: not-allowed;
        }
        .file-info {
            margin-top: 20px;
            padding: 15px;
            background: #f0f2ff;
            border-radius: 6px;
            display: none;
        }
        .file-info.show {
            display: block;
        }
        .file-name {
            color: #333;
            font-weight: bold;
            margin-bottom: 5px;
        }
        .file-size {
            color: #666;
            font-size: 12px;
        }
        .loading-spinner {
            display: none;
            margin-top: 20px;
        }
        .loading-spinner.show {
            display: block;
        }
        .spinner {
            border: 3px solid #f3f3f3;
            border-top: 3px solid #667eea;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto;
        }
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
        .error-message {
            color: #e74c3c;
            margin-top: 15px;
            padding: 10px;
            background: #ffeaea;
            border-radius: 6px;
            display: none;
        }
        .error-message.show {
            display: block;
        }
    </style>
    <!-- Font Awesome for icons -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <!-- jQuery (required for Spectrum) -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <!-- Spectrum Color Picker -->
    <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/spectrum-colorpicker2/dist/spectrum.min.css">
    <script src="https://cdn.jsdelivr.net/npm/spectrum-colorpicker2/dist/spectrum.min.js"></script>
    <style>
        * {
            box-sizing: border-box;
        }
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            background-color: #ffffff;
            overflow: hidden;
        }
        .main-layout {
            position: relative;
            width: 100vw;
            height: 100vh;
            overflow: visible;
            padding-top: 50px;
            padding-bottom: 40px;
        }
        /* Header */
        .main-header {
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 50px;
            background: rgba(44, 62, 80, 0.98);
            backdrop-filter: blur(10px);
            border-bottom: 1px solid rgba(0,0,0,0.1);
            z-index: 3000;
            display: flex;
            align-items: center;
            padding: 0 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .header-brand {
            display: flex;
            align-items: center;
            gap: 10px;
            color: white;
            font-size: 18px;
            font-weight: bold;
            text-decoration: none;
            margin-right: 30px;
        }
        .header-brand i {
            font-size: 24px;
            color: #4a6fa5;
        }
        .header-nav {
            display: flex;
            align-items: center;
            gap: 20px;
            flex: 1;
        }
        .header-nav-item {
            color: rgba(255,255,255,0.9);
            text-decoration: none;
            font-size: 14px;
            padding: 5px 10px;
            border-radius: 4px;
            transition: background 0.2s;
        }
        .header-nav-item:hover {
            background: rgba(255,255,255,0.1);
            color: white;
        }
        .header-right {
            display: flex;
            align-items: center;
            gap: 15px;
            color: rgba(255,255,255,0.8);
            font-size: 12px;
        }
        .header-version {
            padding: 4px 8px;
            background: rgba(255,255,255,0.1);
            border-radius: 4px;
        }
        /* Footer */
        .main-footer {
            position: fixed;
            bottom: 0;
            left: 0;
            width: 100%;
            height: 40px;
            background: rgba(44, 62, 80, 0.98);
            backdrop-filter: blur(10px);
            border-top: 1px solid rgba(0,0,0,0.1);
            z-index: 3000;
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 0 20px;
            box-shadow: 0 -2px 10px rgba(0,0,0,0.1);
            font-size: 12px;
            color: rgba(255,255,255,0.8);
        }
        .footer-left {
            display: flex;
            align-items: center;
            gap: 15px;
        }
        .footer-right {
            display: flex;
            align-items: center;
            gap: 15px;
        }
        .footer-link {
            color: rgba(255,255,255,0.8);
            text-decoration: none;
            transition: color 0.2s;
        }
        .footer-link:hover {
            color: white;
        }
        .footer-separator {
            color: rgba(255,255,255,0.4);
        }
        /* Left sidebar for zoom controls - floating */
        .left-sidebar {
            position: fixed;
            left: 0;
            top: 50px;
            width: 60px;
            height: calc(100vh - 90px);
            background: rgba(44, 62, 80, 0.95);
            backdrop-filter: blur(10px);
            display: flex;
            flex-direction: column;
            align-items: center;
            padding: 15px 5px;
            box-shadow: 2px 0 10px rgba(0,0,0,0.15);
            z-index: 2000;
            border-right: none;
        }
        .zoom-control {
            display: flex;
            flex-direction: column;
            gap: 10px;
            margin-bottom: 20px;
        }
        .zoom-btn {
            width: 40px;
            height: 40px;
            background: #34495e;
            border: none;
            border-radius: 4px;
            color: white;
            font-size: 18px;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            transition: background 0.2s;
            position: relative;
            padding: 0;
            margin: 0 auto;
        }
        .zoom-btn i {
            display: flex;
            align-items: center;
            justify-content: center;
            width: 100%;
            height: 100%;
            margin: 0;
            padding: 0;
        }
        .zoom-btn:hover {
            background: #4a6fa5;
        }
        .zoom-btn:active {
            background: #2c3e50;
        }
        .zoom-value {
            color: white;
            font-size: 11px;
            text-align: center;
            margin-top: 5px;
        }
        .sidebar-divider {
            width: 30px;
            height: 1px;
            background: #34495e;
            margin: 10px 0;
        }
        /* Main content area - full screen */
        .main-content {
            position: absolute;
            top: 50px;
            left: 0;
            width: 100vw;
            height: calc(100vh - 90px);
            overflow: visible;
            background-color: #ffffff;
        }
        .header {
            display: none; /* Hide header for full-screen experience */
        }
        .header h1 {
            margin: 0;
            font-size: 20px;
            color: #333;
        }
        .svg-container {
            position: absolute;
            top: 0;
            left: 0;
            width: 100vw;
            height: calc(100vh - 90px);
            overflow: visible;
            background: #ffffff;
            padding: 0;
            margin: 0;
        }
        .svg-container svg {
            width: 100%;
            height: 100%;
            display: block;
        }
        /* Right floating control panel */
        .control-panel {
            position: fixed;
            top: 70px;
            right: 20px;
            width: 380px;
            max-height: calc(100vh - 110px);
            background: rgba(255, 255, 255, 0.98);
            backdrop-filter: blur(10px);
            border-radius: 8px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.2);
            z-index: 2000;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            user-select: none;
            border: 1px solid rgba(0,0,0,0.1);
        }
        .panel-header {
            background: #007bff;
            color: white;
            padding: 12px 15px;
            font-weight: bold;
            font-size: 14px;
            border-radius: 8px 8px 0 0;
            display: flex;
            justify-content: space-between;
            align-items: center;
            cursor: move;
        }
        .panel-header:active {
            cursor: grabbing;
        }
        .panel-header-buttons {
            display: flex;
            gap: 5px;
            align-items: center;
        }
        .panel-toggle {
            background: rgba(255,255,255,0.2);
            border: none;
            color: white;
            width: 24px;
            height: 24px;
            border-radius: 3px;
            cursor: pointer;
            font-size: 16px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
        }
        .panel-toggle:hover {
            background: rgba(255,255,255,0.3);
        }
        .panel-drag-handle {
            flex: 1;
            cursor: move;
        }
        .panel-content {
            flex: 1;
            overflow-y: auto;
            padding: 0;
        }
        .panel-collapsed .panel-content {
            display: none;
        }
        /* Tabs */
        .tabs {
            display: flex;
            background: #f8f9fa;
            border-bottom: 1px solid #ddd;
        }
        .tab {
            flex: 1;
            padding: 10px 15px;
            background: #f8f9fa;
            border: none;
            border-bottom: 2px solid transparent;
            cursor: pointer;
            font-size: 12px;
            font-weight: 500;
            color: #666;
            transition: all 0.2s;
        }
        .tab:hover {
            background: #e9ecef;
        }
        .tab.active {
            background: white;
            color: #007bff;
            border-bottom-color: #007bff;
        }
        .tab-content {
            display: none;
            padding: 15px;
        }
        .tab-content.active {
            display: block;
        }
        /* Control groups */
        .control-group {
            margin-bottom: 15px;
        }
        .control-group label {
            display: block;
            font-size: 11px;
            font-weight: bold;
            color: #555;
            margin-bottom: 5px;
        }
        .control-group input, .control-group select {
            width: 100%;
            padding: 6px 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 13px;
        }
        .control-group input[type="range"] {
            padding: 0;
        }
        .control-group input[type="checkbox"] {
            width: auto;
            margin-right: 5px;
        }
        .control-group button {
            width: 100%;
            padding: 8px 15px;
            background: #007bff;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 13px;
            margin-top: 5px;
        }
        .control-group button:hover {
            background: #0056b3;
        }
        .control-group button.secondary {
            background: #6c757d;
        }
        .control-group button.secondary:hover {
            background: #5a6268;
        }
        .value-display {
            font-size: 11px;
            color: #333;
            font-weight: bold;
            margin-top: 3px;
        }
        .minimap {
            position: fixed;
            bottom: 60px;
            right: 400px;
            width: 200px;
            height: 150px;
            border: 2px solid #007bff;
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            border-radius: 4px;
            z-index: 1500;
            pointer-events: all;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .minimap svg {
            width: 100%;
            height: 100%;
        }
        .minimap-viewport {
            fill: rgba(0, 123, 255, 0.3);
            stroke: #007bff;
            stroke-width: 2;
            pointer-events: all;
            cursor: move;
        }
        .color-picker-group {
            display: flex;
            flex-direction: column;
            gap: 5px;
            margin-top: 5px;
        }
        .color-picker-item {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .color-picker-item label {
            font-size: 10px;
            width: 80px;
        }
        .color-picker-item input[type="color"] {
            width: 40px;
            height: 25px;
            border: 1px solid #ddd;
            border-radius: 3px;
            cursor: pointer;
        }
        /* Color preview */
        .color-preview {
            display: inline-block;
            width: 30px;
            height: 20px;
            border: 1px solid #ccc;
            border-radius: 3px;
            vertical-align: middle;
            margin-left: 5px;
        }
        /* Context menu */
        .context-menu {
            position: fixed;
            background: white;
            border: 1px solid #ddd;
            border-radius: 4px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            z-index: 10000;
            min-width: 180px;
            display: none;
        }
        .context-menu-header {
            padding: 8px 12px;
            background: #f8f9fa;
            border-bottom: 1px solid #ddd;
            font-weight: bold;
            font-size: 12px;
        }
        .context-menu-item {
            padding: 8px 12px;
            cursor: pointer;
            font-size: 12px;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        .context-menu-item:hover {
            background: #f0f0f0;
        }
        .context-menu-divider {
            height: 1px;
            background: #ddd;
            margin: 4px 0;
        }
        /* Selection box */
        .selection-box {
            position: fixed;
            border: 2px dashed #007bff;
            background: rgba(0, 123, 255, 0.1);
            pointer-events: none;
            z-index: 9999;
        }
        /* Disabled button */
        .zoom-btn:disabled {
            opacity: 0.5;
            cursor: not-allowed;
        }
        .tool-btn.active {
            background: #4a6fa5 !important;
        }
        .tool-selector {
            width: 100%;
        }
        .tooltip {
            position: fixed;
            background: rgba(0, 0, 0, 0.9);
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 12px;
            pointer-events: none;
            z-index: 3000;
            display: none;
            max-width: 250px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.3);
        }
        .tooltip.visible {
            display: block;
        }
        .legend {
            position: fixed;
            bottom: 60px;
            left: 80px;
            padding: 10px 15px;
            background: rgba(249, 249, 249, 0.95);
            backdrop-filter: blur(10px);
            border-radius: 5px;
            font-size: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            z-index: 1500;
            max-width: calc(100vw - 450px);
            border: 1px solid rgba(0,0,0,0.1);
        }
        .legend-item {
            display: inline-block;
            margin-right: 20px;
            margin-bottom: 5px;
        }
        .legend-color {
            display: inline-block;
            width: 20px;
            height: 20px;
            border-radius: 50%;
            margin-right: 5px;
            vertical-align: middle;
            border: 1px solid #ccc;
        }
        .info-panel {
            margin-top: 15px;
            padding: 10px;
            background: #e8f4f8;
            border-left: 4px solid #007bff;
            border-radius: 4px;
            font-size: 13px;
        }
        .highlighted {
            stroke: #ff0000 !important;
            stroke-width: 3px !important;
            opacity: 1 !important;
        }
        .filtered-out {
            opacity: 0.2 !important;
            pointer-events: none;
        }
        .base-control-panel {
            padding: 10px;
            background: #f8f9fa;
            border: 1px solid #ddd;
            border-radius: 4px;
            margin-bottom: 10px;
        }
        .base-header {
            margin-bottom: 8px;
            font-weight: bold;
            font-size: 12px;
        }
        .base-header label {
            display: flex;
            align-items: center;
            gap: 5px;
            cursor: pointer;
        }
        .base-header input[type="checkbox"] {
            width: auto;
            margin: 0;
        }
        .color-range-controls {
            display: flex;
            flex-direction: column;
            gap: 8px;
        }
        .range-control {
            display: flex;
            flex-direction: column;
            gap: 3px;
        }
        .range-control label {
            font-size: 10px;
            color: #666;
        }
        .range-control input[type="range"] {
            width: 100%;
        }
        .range-control span {
            font-size: 10px;
            color: #333;
            font-weight: bold;
        }
        /* Reactivity distribution preview */
        .reactivity-distribution-preview {
            width: 100%;
            height: 60px;
            margin-bottom: 8px;
            border: 1px solid #ddd;
            border-radius: 4px;
            background: transparent;
            position: relative;
            overflow: hidden;
        }
        .reactivity-distribution-preview canvas {
            width: 100%;
            height: 100%;
            display: block;
            background: white;
        }
        .distribution-label {
            font-size: 9px;
            color: #666;
            margin-bottom: 3px;
        }
        .base-disabled {
            opacity: 0.5;
            pointer-events: none;
        }
        .checkbox-group {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .checkbox-group label {
            margin: 0;
            font-weight: normal;
        }
        /* Collapsible sections */
        .control-section {
            margin-bottom: 10px;
        }
        .section-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 8px 10px;
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            cursor: pointer;
            font-weight: 600;
            font-size: 12px;
            color: #495057;
        }
        .section-header:hover {
            background: #e9ecef;
        }
        .section-header i {
            transition: transform 0.2s;
        }
        .section-header.collapsed i {
            transform: rotate(-90deg);
        }
        .section-content {
            padding: 10px;
            border: 1px solid #dee2e6;
            border-top: none;
            border-radius: 0 0 4px 4px;
            background: white;
        }
        .section-content.collapsed {
            display: none;
        }
        /* Table layout for controls */
        .control-table {
            width: 100%;
            border-collapse: collapse;
        }
        .control-table th {
            text-align: left;
            padding: 8px 10px;
            font-size: 11px;
            font-weight: 600;
            color: #495057;
            background: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
        }
        .control-table td {
            padding: 8px 10px;
            border-bottom: 1px solid #f0f0f0;
        }
        .control-table tr:last-child td {
            border-bottom: none;
        }
        /* Tooltip */
        [data-tooltip] {
            position: relative;
        }
        [data-tooltip]:hover::after {
            content: attr(data-tooltip);
            position: absolute;
            left: 100%;
            top: 50%;
            transform: translateY(-50%);
            margin-left: 10px;
            padding: 6px 10px;
            background: rgba(0, 0, 0, 0.9);
            color: white;
            font-size: 11px;
            white-space: nowrap;
            border-radius: 4px;
            z-index: 1000;
            pointer-events: none;
        }
        [data-tooltip]:hover::before {
            content: '';
            position: absolute;
            left: 100%;
            top: 50%;
            transform: translateY(-50%);
            margin-left: 5px;
            border: 5px solid transparent;
            border-right-color: rgba(0, 0, 0, 0.9);
            z-index: 1000;
            pointer-events: none;
        }
        /* Icon buttons */
        .icon-btn {
            background: #007bff;
            color: white;
            border: none;
            border-radius: 4px;
            padding: 8px 12px;
            cursor: pointer;
            display: inline-flex;
            align-items: center;
            gap: 6px;
            font-size: 12px;
            transition: all 0.2s;
        }
        .icon-btn:hover {
            background: #0056b3;
            transform: translateY(-1px);
        }
        .icon-btn i {
            font-size: 14px;
        }
        /* Input with suffix */
        .input-suffix {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .input-suffix input {
            flex: 1;
        }
        .input-suffix span {
            font-size: 11px;
            color: #666;
            white-space: nowrap;
        }
        /* Info panel */
        .info-panel {
            position: fixed;
            top: 110px;
            left: 80px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            padding: 15px;
            max-width: 300px;
            z-index: 998;
            display: none;
        }
        .info-panel.visible {
            display: block;
        }
        .info-panel h3 {
            margin: 0 0 10px 0;
            font-size: 14px;
            color: #333;
        }
        .info-item {
            margin-bottom: 8px;
            font-size: 12px;
        }
        .info-item strong {
            color: #495057;
        }
        /* Search panel */
        .search-panel {
            position: fixed;
            top: 110px;
            left: 80px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            padding: 15px;
            width: 300px;
            z-index: 998;
            display: none;
        }
        .search-panel.visible {
            display: block;
        }
        .search-input {
            width: 100%;
            padding: 8px 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 13px;
        }
        .search-results {
            max-height: 300px;
            overflow-y: auto;
            margin-top: 10px;
        }
        .search-result-item {
            padding: 8px;
            cursor: pointer;
            border-bottom: 1px solid #f0f0f0;
        }
        .search-result-item:hover {
            background: #f8f9fa;
        }
    </style>
    </head>
<body>
    <script>
        // SVG filename for loading message (injected from Rust)
        const EMBEDDED_SVG_FILENAME = "###);
    html.push_str(&format!("'{}';", svg_filename));
    html.push_str(r###"
    </script>
    <!-- Header -->
    <header class="main-header">
        <a href="https://github.com/TongZhou2017/modtector" class="header-brand" target="_blank">
            <i class="fas fa-dna"></i>
            <span>ModTector</span>
        </a>
        <nav class="header-nav">
            <a href="https://modtector.readthedocs.io/" class="header-nav-item" target="_blank" title="Documentation">
                <i class="fas fa-book"></i> Documentation
            </a>
            <a href="https://github.com/TongZhou2017/modtector" class="header-nav-item" target="_blank" title="GitHub Repository">
                <i class="fab fa-github"></i> GitHub
            </a>
            <a href="https://github.com/TongZhou2017/modtector/issues" class="header-nav-item" target="_blank" title="Report Issues">
                <i class="fas fa-bug"></i> Issues
            </a>
        </nav>
        <div class="header-right">
            <span class="header-version">v"###);
    html.push_str(version);
    html.push_str(r###"</span>
        </div>
    </header>
    
    <!-- Initial Page -->
    <div id="initialPage">
        <div class="initial-container">
            <h1><i class="fas fa-dna"></i> RNA Structure Visualizer</h1>
            <p>Please upload an SVG file to start interactive visualization</p>
            
            <div class="file-upload-area" id="uploadArea">
                <div class="file-upload-icon">
                    <i class="fas fa-cloud-upload-alt"></i>
                </div>
                <div class="file-upload-text">Click or drag SVG file here</div>
                <div class="file-upload-hint">Supports .svg format files</div>
            </div>
            
            <input type="file" id="fileInput" accept=".svg">
            <button class="upload-button" id="uploadButton" onclick="document.getElementById('fileInput').click()">
                <i class="fas fa-folder-open"></i> Select File
            </button>
            
            <div class="file-info" id="fileInfo">
                <div class="file-name" id="fileName"></div>
                <div class="file-size" id="fileSize"></div>
            </div>
            
            <div class="loading-spinner" id="loadingSpinner">
                <div class="spinner"></div>
                <p style="margin-top: 10px; color: #666;">Loading...</p>
            </div>
            
            <div class="error-message" id="errorMessage"></div>
        </div>
    </div>
    
    <!-- Main Visualization (hidden initially) -->
    <div id="mainVisualization" style="display: none;">
    <div class="main-layout">
        <!-- Left sidebar for zoom controls -->
        <div class="left-sidebar">
            <div class="zoom-control">
                <button class="zoom-btn" onclick="zoomIn()" data-tooltip="Zoom In"><i class="fas fa-plus"></i></button>
                <button class="zoom-btn" onclick="zoomOut()" data-tooltip="Zoom Out"><i class="fas fa-minus"></i></button>
                <button class="zoom-btn" onclick="resetZoom()" data-tooltip="Fit to Screen"><i class="fas fa-expand-arrows-alt"></i></button>
                <div class="zoom-value" id="zoomValueDisplay">1.0x</div>
            </div>
            <div class="sidebar-divider"></div>
            <button class="zoom-btn" onclick="undoAction()" data-tooltip="Undo (Ctrl+Z)" id="undoBtn" disabled><i class="fas fa-undo"></i></button>
            <button class="zoom-btn" onclick="redoAction()" data-tooltip="Redo (Ctrl+Y)" id="redoBtn" disabled><i class="fas fa-redo"></i></button>
            <div class="sidebar-divider"></div>
            <div class="tool-selector" style="display: flex; flex-direction: column; gap: 5px; margin-bottom: 10px;">
                <button class="zoom-btn tool-btn active" onclick="setToolMode('select')" data-tooltip="Select Tool (Click to select)" id="toolSelect" data-tool="select"><i class="fas fa-mouse-pointer"></i></button>
                <button class="zoom-btn tool-btn" onclick="setToolMode('box')" data-tooltip="Box Select Tool (Drag to select area)" id="toolBox" data-tool="box"><i class="fas fa-vector-square"></i></button>
                <button class="zoom-btn tool-btn" onclick="setToolMode('hand')" data-tooltip="Hand Tool (Pan view)" id="toolHand" data-tool="hand"><i class="fas fa-hand-paper"></i></button>
            </div>
            <div class="sidebar-divider"></div>
            <button class="zoom-btn" onclick="showUploadDialog()" data-tooltip="Upload/Replace SVG"><i class="fas fa-upload"></i></button>
            <div class="sidebar-divider"></div>
            <button class="zoom-btn" onclick="showInfoPanel()" data-tooltip="Information"><i class="fas fa-info-circle"></i></button>
            <button class="zoom-btn" onclick="showSearchPanel()" data-tooltip="Search"><i class="fas fa-search"></i></button>
            <button class="zoom-btn" onclick="showKeyboardHelp()" data-tooltip="Keyboard Shortcuts"><i class="fas fa-keyboard"></i></button>
        </div>
        
        <!-- Main content area -->
        <div class="main-content">
            <div class="header">
                <h1>Interactive RNA Structure Visualization</h1>
            </div>
            <div class="svg-container" id="svgContainer">"###);
    
    // Embed the SVG content  
    html.push_str(svg_content);
    html.push_str(r###"
            <div class="minimap" id="minimap" style="display: none;">
                <svg id="minimapSvg" viewBox="0 0 100 100" preserveAspectRatio="none">
                    <rect class="minimap-viewport" id="minimapViewport" x="0" y="0" width="100" height="100"></rect>
                </svg>
            </div>
        </div>
        <div class="legend" id="legend"></div>
        </div>
        
        <!-- Right floating control panel -->
        <div class="control-panel" id="controlPanel">
            <div class="panel-header" id="panelHeader">
                <span class="panel-drag-handle">Control Panel</span>
                <div class="panel-header-buttons">
                    <button class="panel-toggle" onclick="togglePanel()" id="panelToggle">−</button>
                </div>
            </div>
            <div class="panel-content" id="panelContent">
                <div class="tabs">
                    <button class="tab active" onclick="switchTab('basic')">Basic</button>
                    <button class="tab" onclick="switchTab('threshold')">Threshold</button>
                    <button class="tab" onclick="switchTab('export')">Export</button>
                </div>
                
                <!-- Basic Tab (renamed from Visualization) -->
                <div id="tab-basic" class="tab-content active">
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('background-section')">
                            <span><i class="fas fa-image"></i> Background</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="background-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Background Color</th>
                                    <td>
                                        <div style="display: flex; align-items: center; gap: 5px;">
                                            <input type="text" id="backgroundColor" value="#ffffff" class="spectrum-color-picker" style="flex: 1; padding: 6px 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 13px;">
                                            <span class="color-preview" id="backgroundColorPreview" style="background-color: #ffffff; cursor: pointer; width: 40px; height: 30px; border: 1px solid #ccc; border-radius: 4px; display: inline-block; flex-shrink: 0;"></span>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                    
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('circle-section')">
                            <span><i class="fas fa-circle"></i> Circle Style</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="circle-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Style</th>
                                    <td>
                                        <select id="circleStyle">
                                            <option value="filled" selected>Filled Circles</option>
                                            <option value="hollow">Hollow Rings</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Stroke Width</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="strokeWidth" min="0.5" max="5" step="0.1" value="2" style="flex: 1;">
                                            <span id="strokeWidthValue">2.0</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Circle Radius</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="circleRadius" min="1" max="10" step="0.1" value="2.8" style="flex: 1;">
                                            <span id="circleRadiusValue">2.8</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                    
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('font-section')">
                            <span><i class="fas fa-font"></i> Font Options</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="font-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Font Size</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="fontSize" min="1" max="50" step="0.5" value="12" style="flex: 1;">
                                            <span id="fontSizeValue">12</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Font Color</th>
                                    <td>
                                        <div style="display: flex; align-items: center; gap: 5px;">
                                            <input type="text" id="fontColor" value="#000000" class="spectrum-color-picker" style="flex: 1; padding: 6px 10px; border: 1px solid #ddd; border-radius: 4px; font-size: 13px;">
                                            <span class="color-preview" id="fontColorPreview" style="background-color: #000000; cursor: pointer; width: 40px; height: 30px; border: 1px solid #ccc; border-radius: 4px; display: inline-block; flex-shrink: 0;"></span>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                    
                    <div class="control-group">
                        <button class="icon-btn" onclick="resetView()">
                            <i class="fas fa-redo"></i>
                            <span>Reset View</span>
                        </button>
                    </div>
                </div>
                
                <!-- Threshold Tab (renamed from Bases, merged with Legend) -->
                <div id="tab-threshold" class="tab-content">
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('unified-section')">
                            <span><i class="fas fa-link"></i> Unified Color Ranges</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="unified-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Enable Unified</th>
                                    <td>
                                        <div class="checkbox-group">
                                            <input type="checkbox" id="unifiedColorRanges" checked onchange="toggleUnifiedMode()">
                                            <label for="unifiedColorRanges">Use unified color ranges for all bases</label>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                    
                    <div id="unified-controls" class="control-section" style="display: none;">
                        <div class="section-header" onclick="toggleSection('unified-threshold-section')">
                            <span><i class="fas fa-sliders-h"></i> Unified Thresholds</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="unified-threshold-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Distribution</th>
                                    <td>
                                        <div class="distribution-label">Reactivity Distribution (All Bases)</div>
                                        <div class="reactivity-distribution-preview" id="unifiedDistributionPreview">
                                            <canvas id="unifiedDistributionCanvas"></canvas>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Low-Mid Threshold</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="unifiedRange1" min="0" max="1" step="0.01" value="{unified1_val}" style="flex: 1;">
                                            <span id="unifiedRange1Value">{unified1_str}</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Mid-High Threshold</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="unifiedRange2" min="0" max="1" step="0.01" value="{unified2_val}" style="flex: 1;">
                                            <span id="unifiedRange2Value">{unified2_str}</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Low Color</th>
                                    <td>
                                        <div style="display: flex; align-items: center; gap: 5px;">
                                            <input type="text" id="unifiedColor1" value="#f5f5f5" class="spectrum-color-picker" style="flex: 1;">
                                            <span class="color-preview" id="unifiedColor1Preview" style="background-color: #f5f5f5;"></span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Mid Color</th>
                                    <td>
                                        <div style="display: flex; align-items: center; gap: 5px;">
                                            <input type="text" id="unifiedColor2" value="#ffa500" class="spectrum-color-picker" style="flex: 1;">
                                            <span class="color-preview" id="unifiedColor2Preview" style="background-color: #ffa500;"></span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>High Color</th>
                                    <td>
                                        <div style="display: flex; align-items: center; gap: 5px;">
                                            <input type="text" id="unifiedColor3" value="#ff0000" class="spectrum-color-picker" style="flex: 1;">
                                            <span class="color-preview" id="unifiedColor3Preview" style="background-color: #ff0000;"></span>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                    
                    <div id="individual-controls" class="control-section">"###);
    
    // Set initial values from color_ranges
    let a1_val = color_ranges.a.get(1).map(|r| r.min).unwrap_or(0.3);
    let a2_val = color_ranges.a.get(2).map(|r| r.min).unwrap_or(0.65);
    let t1_val = color_ranges.t.get(1).map(|r| r.min).unwrap_or(0.2);
    let t2_val = color_ranges.t.get(2).map(|r| r.min).unwrap_or(0.55);
    let c1_val = color_ranges.c.get(1).map(|r| r.min).unwrap_or(0.3);
    let c2_val = color_ranges.c.get(2).map(|r| r.min).unwrap_or(0.6);
    let g1_val = color_ranges.g.get(1).map(|r| r.min).unwrap_or(0.2);
    let g2_val = color_ranges.g.get(2).map(|r| r.min).unwrap_or(0.5);
    
    // Unified mode initial values (use A values as default)
    let unified1_val = a1_val;
    let unified2_val = a2_val;
    let unified1_str = format!("{:.2}", unified1_val);
    let unified2_str = format!("{:.2}", unified2_val);
    
    // Replace placeholders in unified controls HTML
    html = html.replace("{unified1_val}", &unified1_val.to_string());
    html = html.replace("{unified1_str}", &unified1_str);
    html = html.replace("{unified2_val}", &unified2_val.to_string());
    html = html.replace("{unified2_str}", &unified2_str);
    
    // Generate base control panels using helper function
    html.push_str(&generate_base_control_panel('A', "Adenine", a1_val, a2_val));
    html.push_str(&generate_base_control_panel('T', "Thymine", t1_val, t2_val));
    html.push_str(&generate_base_control_panel('C', "Cytosine", c1_val, c2_val));
    html.push_str(&generate_base_control_panel('G', "Guanine", g1_val, g2_val));
    
    html.push_str(r###"
                    </div>
                    
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('legend-section')">
                            <span><i class="fas fa-list"></i> Legend Options</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="legend-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>Position</th>
                                    <td>
                                        <select id="legendPosition">
                                            <option value="top">Top</option>
                                            <option value="bottom" selected>Bottom</option>
                                            <option value="left">Left</option>
                                            <option value="right">Right</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Direction</th>
                                    <td>
                                        <select id="legendDirection">
                                            <option value="horizontal" selected>Horizontal</option>
                                            <option value="vertical">Vertical</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Font Size</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendFontSize" min="8" max="24" step="0.5" value="12" style="flex: 1;">
                                            <span id="legendFontSizeValue">12</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Text Alignment</th>
                                    <td>
                                        <select id="legendTextAlign">
                                            <option value="start">Start (Left/Top)</option>
                                            <option value="middle" selected>Middle (Center)</option>
                                            <option value="end">End (Right/Bottom)</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Text Offset X</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendTextOffsetX" min="-50" max="50" step="1" value="0" style="flex: 1;">
                                            <span id="legendTextOffsetXValue">0</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Text Offset Y</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendTextOffsetY" min="-50" max="50" step="1" value="0" style="flex: 1;">
                                            <span id="legendTextOffsetYValue">0</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Base Spacing (Rows)</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendRowSpacing" min="0" max="100" step="1" value="40" style="flex: 1;">
                                            <span id="legendRowSpacingValue">40</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Base Spacing (Columns)</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendColSpacing" min="0" max="100" step="1" value="50" style="flex: 1;">
                                            <span id="legendColSpacingValue">50</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Item Width</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendItemWidth" min="10" max="100" step="1" value="30" style="flex: 1;">
                                            <span id="legendItemWidthValue">30</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                                <tr>
                                    <th>Item Height</th>
                                    <td>
                                        <div class="input-suffix">
                                            <input type="range" id="legendItemHeight" min="5" max="50" step="1" value="15" style="flex: 1;">
                                            <span id="legendItemHeightValue">15</span>
                                            <span>px</span>
                                        </div>
                                    </td>
                                </tr>
                            </table>
                        </div>
                    </div>
                </div>
                
                <!-- Export Tab -->
                <div id="tab-export" class="tab-content">
                    <div class="control-section">
                        <div class="section-header" onclick="toggleSection('export-section')">
                            <span><i class="fas fa-download"></i> Export Options</span>
                            <i class="fas fa-chevron-down"></i>
                        </div>
                        <div id="export-section" class="section-content">
                            <table class="control-table">
                                <tr>
                                    <th>File Name</th>
                                    <td>
                                        <input type="text" id="exportFileName" value="rna_structure_interactive" style="width: 100%; padding: 6px;">
                                    </td>
                                </tr>
                                <tr>
                                    <th>Format</th>
                                    <td>
                                        <select id="exportFormat" style="width: 100%; padding: 6px;">
                                            <option value="svg" selected>SVG</option>
                                            <option value="png">PNG</option>
                                            <option value="pdf">PDF</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr id="pngResolutionRow" style="display: none;">
                                    <th>PNG Resolution</th>
                                    <td>
                                        <select id="pngResolution" style="width: 100%; padding: 6px;">
                                            <option value="1" selected>1x (Standard)</option>
                                            <option value="2">2x (High)</option>
                                            <option value="3">3x (Very High)</option>
                                        </select>
                                    </td>
                                </tr>
                                <tr id="pngDPIRow" style="display: none;">
                                    <th>DPI</th>
                                    <td>
                                        <select id="pngDPI" style="width: 100%; padding: 6px;">
                                            <option value="72" selected>72 DPI</option>
                                            <option value="150">150 DPI</option>
                                            <option value="300">300 DPI</option>
                                            <option value="600">600 DPI</option>
                                        </select>
                                    </td>
                                </tr>
                            </table>
                            <div style="margin-top: 15px;">
                                <button class="icon-btn" onclick="exportImage()" style="width: 100%;">
                                    <i class="fas fa-download"></i>
                                    <span>Export</span>
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    </div>
    
    <div class="tooltip" id="tooltip"></div>
    
    <!-- Footer -->
    <footer class="main-footer">
        <div class="footer-left">
            <span><strong>ModTector</strong> - A high-performance RNA modification detection tool</span>
            <span class="footer-separator">|</span>
            <span>Version "###);
    html.push_str(version);
    html.push_str(r###"</span>
            <span class="footer-separator">|</span>
            <a href="https://github.com/TongZhou2017/modtector" class="footer-link" target="_blank">GitHub</a>
            <span class="footer-separator">|</span>
            <a href="https://modtector.readthedocs.io/" class="footer-link" target="_blank">Documentation</a>
        </div>
        <div class="footer-right">
            <span>License: MIT</span>
            <span class="footer-separator">|</span>
            <span>Author: Tong Zhou</span>
            <span class="footer-separator">|</span>
            <a href="https://github.com/TongZhou2017/modtector/issues" class="footer-link" target="_blank">Report Issues</a>
        </div>
    </footer>
    
    <!-- Context menu for individual base editing -->
    <div id="contextMenu" class="context-menu" style="display: none;">
        <div class="context-menu-header">Edit Base</div>
        <div class="context-menu-item" onclick="editBaseFontSize()">
            <i class="fas fa-font"></i> Font Size
        </div>
        <div class="context-menu-item" onclick="editBaseFontColor()">
            <i class="fas fa-palette"></i> Font Color
        </div>
        <div class="context-menu-item" onclick="editBaseCircleRadius()">
            <i class="fas fa-circle"></i> Circle Radius
        </div>
        <div class="context-menu-item" onclick="editBaseCircleColor()">
            <i class="fas fa-fill"></i> Circle Color
        </div>
        <div class="context-menu-divider"></div>
        <div class="context-menu-item" onclick="resetBaseStyle()">
            <i class="fas fa-undo"></i> Reset to Default
        </div>
    </div>
    
    <!-- Batch edit menu for multiple bases -->
    <div id="batchEditMenu" class="context-menu" style="display: none;">
        <div class="context-menu-header">Edit Selected Bases (<span id="batchEditCount">0</span>)</div>
        <div class="context-menu-item" onclick="batchEditFontSize()">
            <i class="fas fa-font"></i> Font Size
        </div>
        <div class="context-menu-item" onclick="batchEditFontColor()">
            <i class="fas fa-palette"></i> Font Color
        </div>
        <div class="context-menu-item" onclick="batchEditCircleRadius()">
            <i class="fas fa-circle"></i> Circle Radius
        </div>
        <div class="context-menu-item" onclick="batchEditCircleColor()">
            <i class="fas fa-fill"></i> Circle Color
        </div>
        <div class="context-menu-divider"></div>
        <div class="context-menu-item" onclick="batchResetStyle()">
            <i class="fas fa-undo"></i> Reset to Default
        </div>
    </div>
    
    <!-- Selection box for batch editing -->
    <div id="selectionBox" class="selection-box" style="display: none;"></div>
    
    <!-- Info Panel -->
    <div class="info-panel" id="infoPanel">
        <h3><i class="fas fa-info-circle"></i> Information</h3>
        <div class="info-item">
            <strong>Total Positions:</strong> <span id="infoTotalPositions">0</span>
        </div>
        <div class="info-item">
            <strong>Signal Type:</strong> <span id="infoSignalType">-</span>
        </div>
        <div class="info-item">
            <strong>Visible Bases:</strong> <span id="infoVisibleBases">-</span>
        </div>
        <div class="info-item">
            <strong>Zoom Level:</strong> <span id="infoZoomLevel">1.0x</span>
        </div>
    </div>
    
    <!-- Search Panel -->
    <div class="search-panel" id="searchPanel">
        <h3><i class="fas fa-search"></i> Search</h3>
        <input type="text" class="search-input" id="searchInput" placeholder="Search by position or base..." onkeyup="performSearch()">
        <div class="search-results" id="searchResults"></div>
    </div>
    
    <script>
        // File upload and initialization
        let uploadedSvgContent = null;
        let uploadedReactivityData = null;
        let svgReplaceInput = null;
        
        // Create hidden file input for sidebar button
        function createSvgReplaceInput() {
            if (!svgReplaceInput) {
                svgReplaceInput = document.createElement('input');
                svgReplaceInput.type = 'file';
                svgReplaceInput.accept = '.svg';
                svgReplaceInput.style.display = 'none';
                svgReplaceInput.addEventListener('change', function(e) {
                    if (e.target.files.length > 0) {
                        handleFile(e.target.files[0], true);
                    }
                });
                document.body.appendChild(svgReplaceInput);
            }
            return svgReplaceInput;
        }
        
        // Show upload dialog (for sidebar button)
        function showUploadDialog() {
            const input = createSvgReplaceInput();
            input.click();
        }
        
        // Check if SVG content is already embedded (from Rust generation)
        const existingSvg = document.querySelector('#svgContainer svg');
        if (existingSvg && existingSvg.innerHTML.trim().length > 0) {
            // SVG already loaded, show loading message first
            const initialPage = document.getElementById('initialPage');
            const initialContainer = initialPage.querySelector('.initial-container');
            const uploadArea = document.getElementById('uploadArea');
            const loadingSpinner = document.getElementById('loadingSpinner');
            
            // Get SVG filename from embedded constant or extract from URL
            let svgFileName = typeof EMBEDDED_SVG_FILENAME !== 'undefined' ? EMBEDDED_SVG_FILENAME : 'RNA structure.svg';
            if (svgFileName === 'RNA structure.svg') {
                try {
                    const currentPath = window.location.pathname;
                    const pathParts = currentPath.split('/');
                    const fileName = pathParts[pathParts.length - 1];
                    if (fileName && fileName.endsWith('.html')) {
                        // Try to extract SVG name from HTML filename
                        // e.g., "br_k_m2018_mut_shapemap_EC_1M7_rep2_EC_16S.html" -> "EC_16S.svg"
                        const geneMatch = fileName.match(/(EC_\d+S|Human_\d+S|[\w_]+)/);
                        if (geneMatch) {
                            svgFileName = geneMatch[1] + '.svg';
                        }
                    }
                } catch (e) {
                    console.warn('Could not extract SVG filename:', e);
                }
            }
            
            // Update initial page to show loading message
            if (uploadArea) {
                uploadArea.style.display = 'none';
            }
            if (initialContainer) {
                const loadingText = initialContainer.querySelector('p');
                if (loadingText) {
                    loadingText.textContent = `Loading ${svgFileName}...`;
                    loadingText.style.fontSize = '18px';
                    loadingText.style.color = '#333';
                    loadingText.style.marginBottom = '20px';
                }
            }
            if (loadingSpinner) {
                loadingSpinner.classList.add('show');
                const spinnerText = loadingSpinner.querySelector('p');
                if (spinnerText) {
                    spinnerText.textContent = `Loading ${svgFileName}...`;
                }
            }
            
            // Show loading message for 2.5 seconds before hiding initial page
            setTimeout(() => {
                initialPage.classList.add('hidden');
                setTimeout(() => {
                    initialPage.style.display = 'none';
                    document.getElementById('mainVisualization').style.display = 'block';
                }, 500);
            }, 2500);
        } else {
            // No SVG content, show initial page
            document.getElementById('mainVisualization').style.display = 'none';
        }
        
        // File upload handlers
        const uploadArea = document.getElementById('uploadArea');
        const fileInput = document.getElementById('fileInput');
        const fileInfo = document.getElementById('fileInfo');
        const fileName = document.getElementById('fileName');
        const fileSize = document.getElementById('fileSize');
        const loadingSpinner = document.getElementById('loadingSpinner');
        const errorMessage = document.getElementById('errorMessage');
        
        // Click to upload
        uploadArea.addEventListener('click', () => {
            fileInput.click();
        });
        
        // Drag and drop
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('dragover');
        });
        
        uploadArea.addEventListener('dragleave', () => {
            uploadArea.classList.remove('dragover');
        });
        
        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('dragover');
            const files = e.dataTransfer.files;
            if (files.length > 0) {
                handleFile(files[0]);
            }
        });
        
        // File input change
        fileInput.addEventListener('change', (e) => {
            if (e.target.files.length > 0) {
                handleFile(e.target.files[0]);
            }
        });
        
        // Handle file upload
        function handleFile(file, isReplace = false) {
            // Validate file type
            if (!file.name.toLowerCase().endsWith('.svg')) {
                showError('Please upload an SVG format file');
                return;
            }
            
            // Show file info (only on initial page)
            if (!isReplace) {
                fileName.textContent = file.name;
                fileSize.textContent = formatFileSize(file.fileSize || file.size);
                fileInfo.classList.add('show');
                errorMessage.classList.remove('show');
                
                // Show loading
                loadingSpinner.classList.add('show');
            } else {
                // Show loading indicator for replace operation
                showLoadingIndicator('Loading new SVG file...');
            }
            
            // Read file
            const reader = new FileReader();
            reader.onload = function(e) {
                try {
                    const svgContent = e.target.result;
                    uploadedSvgContent = svgContent;
                    
                    // Parse SVG to check if it's valid
                    const parser = new DOMParser();
                    const svgDoc = parser.parseFromString(svgContent, 'image/svg+xml');
                    const parseError = svgDoc.querySelector('parsererror');
                    
                    if (parseError) {
                        showError('SVG file format error, please check file content');
                        if (!isReplace) {
                            loadingSpinner.classList.remove('show');
                        } else {
                            hideLoadingIndicator();
                        }
                        return;
                    }
                    
                    // Load SVG into container
                    loadSvgIntoVisualization(svgContent, isReplace);
                    
                } catch (error) {
                    showError('Failed to read file: ' + error.message);
                    if (!isReplace) {
                        loadingSpinner.classList.remove('show');
                    } else {
                        hideLoadingIndicator();
                    }
                }
            };
            
            reader.onerror = function() {
                showError('Error occurred while reading file');
                if (!isReplace) {
                    loadingSpinner.classList.remove('show');
                } else {
                    hideLoadingIndicator();
                }
            };
            
            reader.readAsText(file);
        }
        
        // Load SVG into visualization
        function loadSvgIntoVisualization(svgContent, isReplace = false) {
            try {
                // Parse SVG
                const parser = new DOMParser();
                const svgDoc = parser.parseFromString(svgContent, 'image/svg+xml');
                const svgElement = svgDoc.documentElement;
                
                // Get SVG container
                const svgContainer = document.getElementById('svgContainer');
                if (!svgContainer) {
                    showError('SVG container not found');
                    if (!isReplace) {
                        loadingSpinner.classList.remove('show');
                    } else {
                        hideLoadingIndicator();
                    }
                    return;
                }
                
                // Clear existing content and remove old event listeners
                const oldSvg = svgContainer.querySelector('svg');
                if (oldSvg) {
                    // Clone to remove all event listeners
                    const newSvg = oldSvg.cloneNode(false);
                    svgContainer.innerHTML = '';
                } else {
                    svgContainer.innerHTML = '';
                }
                
                // Clone and append SVG
                const clonedSvg = svgElement.cloneNode(true);
                svgContainer.appendChild(clonedSvg);
                
                // Initialize reactivity data if not available
                if (!window.reactivityData || !window.reactivityData.positions) {
                    // Try to extract reactivity data from SVG
                    initializeReactivityDataFromSvg(clonedSvg);
                }
                
                // Initialize viewBox
                const svg = clonedSvg;
                svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
                
                let viewBox = svg.viewBox.baseVal;
                if (viewBox.width === 0 || viewBox.height === 0) {
                    // If viewBox is not set, use SVG dimensions
                    const svgWidth = parseFloat(svg.getAttribute('width')) || svg.getBBox().width;
                    const svgHeight = parseFloat(svg.getAttribute('height')) || svg.getBBox().height;
                    if (svgWidth > 0 && svgHeight > 0) {
                        svg.setAttribute('viewBox', `0 0 ${svgWidth} ${svgHeight}`);
                        viewBox = svg.viewBox.baseVal;
                    }
                }
                
                // Store initial viewBox
                window.initialViewBox = {
                    x: viewBox.x,
                    y: viewBox.y,
                    width: viewBox.width,
                    height: viewBox.height
                };
                
                if (isReplace) {
                    // Replace mode: reinitialize visualization
                    hideLoadingIndicator();
                    
                    // Clear old state
                    selectedBases.clear();
                    highlightedCircles.clear();
                    isSelecting = false;
                    isPanning = false;
                    
                    // Hide any open menus
                    const batchMenu = document.getElementById('batchEditMenu');
                    const contextMenu = document.getElementById('contextMenu');
                    if (batchMenu) batchMenu.style.display = 'none';
                    if (contextMenu) contextMenu.style.display = 'none';
                    
                    // Clear saved state when replacing SVG to use original SVG values
                    // This ensures the new SVG uses its original font sizes and legend positions
                    localStorage.removeItem('modDetectorState');
                    
                    // Reset UI controls to defaults
                    const fontSizeInput = document.getElementById('fontSize');
                    if (fontSizeInput) {
                        fontSizeInput.value = '';
                    }
                    const legendPositionSelect = document.getElementById('legendPosition');
                    if (legendPositionSelect) {
                        legendPositionSelect.value = 'bottom';
                    }
                    
                    // Reinitialize visualization with full setup
                    setTimeout(() => {
                        // Reinitialize SVG interactivity first
                        initializeSvgInteractivity();
                        // Then initialize all other features (will use original SVG values)
                        initializeVisualizationFeatures();
                    }, 100);
                } else {
                    // Initial load: Hide initial page and show visualization
                    setTimeout(() => {
                        document.getElementById('initialPage').classList.add('hidden');
                        setTimeout(() => {
                            document.getElementById('initialPage').style.display = 'none';
                            document.getElementById('mainVisualization').style.display = 'block';
                            
                            // Initialize visualization
                            initializeSvgInteractivity();
                            initializeVisualizationFeatures();
                            
                            loadingSpinner.classList.remove('show');
                        }, 500);
                    }, 500);
                }
                
            } catch (error) {
                showError('Failed to load SVG: ' + error.message);
                console.error('Error loading SVG:', error);
                if (!isReplace) {
                    loadingSpinner.classList.remove('show');
                } else {
                    hideLoadingIndicator();
                }
            }
        }
        
        // Initialize reactivity data from SVG
        function initializeReactivityDataFromSvg(svg) {
            if (!window.reactivityData) {
                window.reactivityData = { signalType: 'unknown', positions: {} };
            }
            
            // Extract data from circles
            const circles = svg.querySelectorAll('circle');
            circles.forEach(circle => {
                const position = circle.getAttribute('data-position');
                if (position) {
                    const base = circle.getAttribute('data-base') || 'N';
                    const reactivity = parseFloat(circle.getAttribute('data-reactivity')) || 0;
                    
                    if (!window.reactivityData.positions[position]) {
                        window.reactivityData.positions[position] = {
                            base: base,
                            reactivity: reactivity
                        };
                    }
                }
            });
        }
        
        // Show error message
        function showError(message) {
            errorMessage.textContent = message;
            errorMessage.classList.add('show');
        }
        
        // Format file size
        function formatFileSize(bytes) {
            if (bytes === 0) return '0 Bytes';
            const k = 1024;
            const sizes = ['Bytes', 'KB', 'MB', 'GB'];
            const i = Math.floor(Math.log(bytes) / Math.log(k));
            return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
        }
        
        // Show loading indicator (for replace operation)
        function showLoadingIndicator(message) {
            let loadingDiv = document.getElementById('replaceLoadingIndicator');
            if (!loadingDiv) {
                loadingDiv = document.createElement('div');
                loadingDiv.id = 'replaceLoadingIndicator';
                loadingDiv.style.cssText = 'position: fixed; top: 20px; right: 20px; background: #667eea; color: white; padding: 15px 25px; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.2); z-index: 10001; display: flex; align-items: center; gap: 10px;';
                document.body.appendChild(loadingDiv);
            }
            loadingDiv.innerHTML = '<div class="spinner" style="width: 20px; height: 20px; border-width: 2px;"></div><span>' + (message || 'Loading...') + '</span>';
            loadingDiv.style.display = 'flex';
        }
        
        // Hide loading indicator
        function hideLoadingIndicator() {
            const loadingDiv = document.getElementById('replaceLoadingIndicator');
            if (loadingDiv) {
                loadingDiv.style.display = 'none';
            }
        }
        
        // Initialize visualization (called after SVG is loaded)
        function initializeVisualization() {
            // This will be called after SVG is loaded
            // Re-run the initialization code
            initializeSvgInteractivity();
            if (typeof initializeVisualizationFeatures === 'function') {
                initializeVisualizationFeatures();
            } else {
                // Fallback: trigger existing initialization
                setTimeout(() => {
                    const event = new Event('svgLoaded');
                    document.dispatchEvent(event);
                }, 100);
            }
        }
        
        // Tab switching
        function switchTab(tabName) {
            // Hide all tabs
            document.querySelectorAll('.tab').forEach(tab => tab.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));
            
            // Show selected tab - map old names to new names for compatibility
            const tabNameMap = {
                'visualization': 'basic',
                'bases': 'threshold',
                'legend': 'threshold'
            };
            const actualTabName = tabNameMap[tabName] || tabName;
            
            // Find and activate the clicked tab button
            const tabs = document.querySelectorAll('.tab');
            tabs.forEach(tab => {
                const onclickStr = tab.getAttribute('onclick') || '';
                if (onclickStr.includes(`'${tabName}'`) || onclickStr.includes(`'${actualTabName}'`)) {
                    tab.classList.add('active');
                }
            });
            
            // Show corresponding content
            const content = document.getElementById('tab-' + actualTabName);
            if (content) {
                content.classList.add('active');
            }
        }
        
        // Panel toggle
        function togglePanel() {
            const panel = document.getElementById('controlPanel');
            const toggle = document.getElementById('panelToggle');
            panel.classList.toggle('panel-collapsed');
            toggle.textContent = panel.classList.contains('panel-collapsed') ? '+' : '−';
        }
        
        // Section toggle (collapse/expand)
        function toggleSection(sectionId) {
            const section = document.getElementById(sectionId);
            const header = section.previousElementSibling;
            section.classList.toggle('collapsed');
            header.classList.toggle('collapsed');
        }
        
        // Info panel
        function showInfoPanel() {
            const panel = document.getElementById('infoPanel');
            const searchPanel = document.getElementById('searchPanel');
            searchPanel.classList.remove('visible');
            panel.classList.toggle('visible');
            updateInfoPanel();
        }
        
        function updateInfoPanel() {
            const totalPositions = Object.keys(reactivityData.positions || {}).length;
            const signalType = reactivityData.signalType || '-';
            const visibleBases = ['A', 'T', 'C', 'G'].filter(base => 
                document.getElementById('baseEnable' + base)?.checked
            ).join(', ');
            const zoomLevel = document.getElementById('zoomLevelValue')?.textContent || '1.0x';
            
            document.getElementById('infoTotalPositions').textContent = totalPositions;
            document.getElementById('infoSignalType').textContent = signalType;
            document.getElementById('infoVisibleBases').textContent = visibleBases || 'None';
            document.getElementById('infoZoomLevel').textContent = zoomLevel;
        }
        
        // Search panel
        function showSearchPanel() {
            const panel = document.getElementById('searchPanel');
            const infoPanel = document.getElementById('infoPanel');
            infoPanel.classList.remove('visible');
            panel.classList.toggle('visible');
            if (panel.classList.contains('visible')) {
                document.getElementById('searchInput').focus();
            }
        }
        
        function performSearch() {
            const query = document.getElementById('searchInput').value.toLowerCase();
            const results = document.getElementById('searchResults');
            results.innerHTML = '';
            
            if (!query) {
                return;
            }
            
            const matches = [];
            for (const [pos, data] of Object.entries(reactivityData.positions || {})) {
                if (pos.includes(query) || data.base.toLowerCase() === query) {
                    matches.push({ position: pos, data: data });
                }
            }
            
            if (matches.length === 0) {
                results.innerHTML = '<div style="padding: 10px; color: #666;">No matches found</div>';
                return;
            }
            
            matches.slice(0, 20).forEach(match => {
                const item = document.createElement('div');
                item.className = 'search-result-item';
                item.innerHTML = `<strong>Position ${match.position}</strong>: ${match.data.base} (${match.data.reactivity.toFixed(4)})`;
                item.onclick = () => {
                    highlightPosition(match.position);
                    showSearchPanel(); // Close panel
                };
                results.appendChild(item);
            });
        }
        
        function highlightPosition(position) {
            const circle = document.querySelector(`circle[data-position="${position}"]`);
            if (circle) {
                circle.scrollIntoView({ behavior: 'smooth', block: 'center' });
                circle.classList.add('highlighted');
                setTimeout(() => circle.classList.remove('highlighted'), 2000);
            }
        }
        
        // Keyboard shortcuts
        function showKeyboardHelp() {
            alert('Keyboard Shortcuts:\n\n' +
                  '+ / = : Zoom In\n' +
                  '- / _ : Zoom Out\n' +
                  '0 : Reset Zoom\n' +
                  'F : Fit to Screen\n' +
                  'R : Reset View\n' +
                  'E : Export SVG\n' +
                  'I : Show Info\n' +
                  'S : Search\n' +
                  'H : Show Help\n' +
                  'Esc : Close Panels');
        }
        
        // Initialize Spectrum color pickers
        function initColorPickers() {
            if (typeof $ === 'undefined') {
                console.warn('jQuery not loaded, color pickers will use default inputs');
                return;
            }
            document.querySelectorAll('.spectrum-color-picker').forEach(input => {
                try {
                    const inputId = input.id;
                    const isBackgroundColor = (inputId === 'backgroundColor');
                    const isUnifiedColor = (inputId === 'unifiedColor1' || inputId === 'unifiedColor2' || inputId === 'unifiedColor3');
                    
                    // Extract base from color picker ID (e.g., 'colorA1' -> 'A')
                    let baseFromId = null;
                    if (!isBackgroundColor && inputId !== 'fontColor' && !isUnifiedColor) {
                        const match = inputId.match(/^color([ATCG])(\d+)$/);
                        if (match) {
                            baseFromId = match[1];
                        }
                    }
                    
                    $(input).spectrum({
                        color: input.value,
                        showInput: true,
                        showAlpha: false,
                        className: "full-spectrum",
                        showInitial: true,
                        showPalette: true,
                        showSelectionPalette: true,
                        maxSelectionSize: 10,
                        preferredFormat: "hex",
                        allowEmpty: false,
                        clickoutFiresChange: true,
                        palette: [
                            ["#000","#444","#666","#999","#ccc","#eee","#f3f3f3","#fff"],
                            ["#f00","#f90","#ff0","#0f0","#0ff","#00f","#90f","#f0f"],
                            ["#f4cccc","#fce5cd","#fff2cc","#d9ead3","#d0e0e3","#cfe2f3","#d9d2e9","#ead1dc"],
                            ["#ea9999","#f9cb9c","#ffe599","#b6d7a8","#a2c4c9","#a4c2f4","#b4a7d6","#d5a6bd"],
                            ["#e06666","#f6b26b","#ffd966","#93c47d","#76a5af","#6d9eeb","#8e7cc3","#c27ba0"],
                            ["#cc0000","#e69138","#f1c232","#6aa84f","#45818e","#3c78d8","#674ea7","#a64d79"],
                            ["#990000","#b45f06","#bf9000","#38761d","#134f5c","#1155cc","#351c75","#741b47"],
                            ["#660000","#783f04","#7f6000","#274e13","#0c343d","#0b5394","#20124d","#4c1130"]
                        ],
                        change: function(color) {
                            if (isBackgroundColor) {
                                updateBackgroundColor();
                                updateColorPreview('backgroundColor', 'backgroundColorPreview');
                                saveState();
                            } else if (inputId === 'fontColor') {
                                updateFontStyle();
                                updateColorPreview('fontColor', 'fontColorPreview');
                                saveState();
                            } else if (isUnifiedColor) {
                                // Unified color picker: update all bases
                                updateColorRanges();
                            } else if (baseFromId) {
                                // CRITICAL FIX: Individual base color picker - only update this specific base
                                updateColorRangeForBase(baseFromId);
                                updateLegendValues(baseFromId); // Only update legend for this base
                                updateColorsForBase(baseFromId); // Only update colors for this base
                                updateDistributionPreviews(baseFromId); // Only update preview for this base
                            } else {
                                // Fallback: update all (should not happen for base-specific colors)
                                updateColorRanges();
                            }
                        }
                    });
                } catch (e) {
                    console.warn('Failed to initialize color picker:', e);
                }
            });
        }
        
        // Load saved state from localStorage
        function loadSavedState() {
            try {
                const saved = localStorage.getItem('rnaVizState');
                if (saved) {
                    const state = JSON.parse(saved);
                    if (state.zoomLevel) document.getElementById('zoomLevel').value = state.zoomLevel;
                    if (state.circleStyle) document.getElementById('circleStyle').value = state.circleStyle;
                    if (state.strokeWidth) document.getElementById('strokeWidth').value = state.strokeWidth;
                    if (state.fontSize) document.getElementById('fontSize').value = state.fontSize;
                    if (state.fontColor) {
                        const fontColorEl = document.getElementById('fontColor');
                        if (fontColorEl) {
                            fontColorEl.value = state.fontColor;
                            if (typeof $ !== 'undefined' && $(fontColorEl).spectrum) {
                                $(fontColorEl).spectrum('set', state.fontColor);
                            }
                        }
                    }
                    if (state.backgroundColor) {
                        const backgroundColorEl = document.getElementById('backgroundColor');
                        if (backgroundColorEl) {
                            backgroundColorEl.value = state.backgroundColor;
                            if (typeof $ !== 'undefined' && $(backgroundColorEl).spectrum) {
                                $(backgroundColorEl).spectrum('set', state.backgroundColor);
                            }
                            updateBackgroundColor();
                        }
                    }
                    if (state.unifiedColorRanges !== undefined) {
                        document.getElementById('unifiedColorRanges').checked = state.unifiedColorRanges;
                    }
                    // Load base enables
                    ['A', 'T', 'C', 'G'].forEach(base => {
                        if (state[`baseEnable${base}`] !== undefined) {
                            document.getElementById(`baseEnable${base}`).checked = state[`baseEnable${base}`];
                        }
                    });
                    // Load color ranges
                    ['A', 'T', 'C', 'G'].forEach(base => {
                        for (let i = 1; i <= 2; i++) {
                            if (state[`range${base}${i}`] !== undefined) {
                                document.getElementById(`range${base}${i}`).value = state[`range${base}${i}`];
                            }
                        }
                        for (let i = 1; i <= 3; i++) {
                            if (state[`color${base}${i}`] !== undefined) {
                                const input = document.getElementById(`color${base}${i}`);
                                if (input) {
                                    input.value = state[`color${base}${i}`];
                                    $(input).spectrum('set', state[`color${base}${i}`]);
                                }
                            }
                        }
                    });
                }
            } catch (e) {
                console.error('Error loading saved state:', e);
            }
        }
        
        // Save state to localStorage
        function saveState() {
            try {
                const state = {};
                const zoomLevelEl = document.getElementById('zoomLevel');
                if (zoomLevelEl) state.zoomLevel = zoomLevelEl.value;
                
                const circleStyleEl = document.getElementById('circleStyle');
                if (circleStyleEl) state.circleStyle = circleStyleEl.value;
                
                const strokeWidthEl = document.getElementById('strokeWidth');
                if (strokeWidthEl) state.strokeWidth = strokeWidthEl.value;
                
                const fontSizeEl = document.getElementById('fontSize');
                if (fontSizeEl) state.fontSize = fontSizeEl.value;
                
                const fontColorEl = document.getElementById('fontColor');
                if (fontColorEl) state.fontColor = fontColorEl.value;
                
                const backgroundColorEl = document.getElementById('backgroundColor');
                if (backgroundColorEl) state.backgroundColor = backgroundColorEl.value;
                
                const unifiedColorRangesEl = document.getElementById('unifiedColorRanges');
                if (unifiedColorRangesEl) state.unifiedColorRanges = unifiedColorRangesEl.checked;
                
                ['A', 'T', 'C', 'G'].forEach(base => {
                    const baseEnableEl = document.getElementById(`baseEnable${base}`);
                    if (baseEnableEl) {
                        state[`baseEnable${base}`] = baseEnableEl.checked;
                    }
                    for (let i = 1; i <= 2; i++) {
                        const rangeEl = document.getElementById(`range${base}${i}`);
                        if (rangeEl) {
                            state[`range${base}${i}`] = rangeEl.value;
                        }
                    }
                    for (let i = 1; i <= 3; i++) {
                        const input = document.getElementById(`color${base}${i}`);
                        if (input) {
                            state[`color${base}${i}`] = input.value;
                        }
                    }
                });
                localStorage.setItem('rnaVizState', JSON.stringify(state));
            } catch (e) {
                console.error('Error saving state:', e);
            }
        }
        
        // Zoom functions for sidebar buttons
        function zoomIn() {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            // Use window.initialViewBox if available, otherwise try to get it
            const initViewBox = window.initialViewBox || initialViewBox;
            if (!initViewBox || initViewBox.width === 0) {
                console.warn('Initial viewBox not set, cannot zoom');
                return;
            }
            
            const viewBox = svg.viewBox.baseVal;
            const rect = svg.getBoundingClientRect();
            if (rect.width === 0 || rect.height === 0) return;
            
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const x = (centerX / rect.width) * viewBox.width + viewBox.x;
            const y = (centerY / rect.height) * viewBox.height + viewBox.y;
            
            const delta = 0.9; // Zoom in
            const newWidth = viewBox.width * delta;
            const newHeight = viewBox.height * delta;
            
            viewBox.x = x - (x - viewBox.x) * delta;
            viewBox.y = y - (y - viewBox.y) * delta;
            viewBox.width = newWidth;
            viewBox.height = newHeight;
            
            // Constrain to bounds
            viewBox.x = Math.max(initViewBox.x, Math.min(initViewBox.x + initViewBox.width - viewBox.width, viewBox.x));
            viewBox.y = Math.max(initViewBox.y, Math.min(initViewBox.y + initViewBox.height - viewBox.height, viewBox.y));
            
            updateMinimap();
            updateZoomDisplayFromViewBox();
            // Reset minimap auto-hide timer on zoom
            resetMinimapHideTimer();
        }
        
        function zoomOut() {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            // Use window.initialViewBox if available, otherwise try to get it
            const initViewBox = window.initialViewBox || initialViewBox;
            if (!initViewBox || initViewBox.width === 0) {
                console.warn('Initial viewBox not set, cannot zoom');
                return;
            }
            
            const viewBox = svg.viewBox.baseVal;
            const rect = svg.getBoundingClientRect();
            if (rect.width === 0 || rect.height === 0) return;
            
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const x = (centerX / rect.width) * viewBox.width + viewBox.x;
            const y = (centerY / rect.height) * viewBox.height + viewBox.y;
            
            const delta = 1.1; // Zoom out
            const newWidth = viewBox.width * delta;
            const newHeight = viewBox.height * delta;
            
            // Don't zoom out beyond initial size
            if (newWidth > initViewBox.width || newHeight > initViewBox.height) {
                resetZoom();
                return;
            }
            
            viewBox.x = x - (x - viewBox.x) * delta;
            viewBox.y = y - (y - viewBox.y) * delta;
            viewBox.width = newWidth;
            viewBox.height = newHeight;
            
            // Constrain to bounds
            viewBox.x = Math.max(initViewBox.x, Math.min(initViewBox.x + initViewBox.width - viewBox.width, viewBox.x));
            viewBox.y = Math.max(initViewBox.y, Math.min(initViewBox.y + initViewBox.height - viewBox.height, viewBox.y));
            
            updateMinimap();
            updateZoomDisplayFromViewBox();
            // Reset minimap auto-hide timer on zoom
            resetMinimapHideTimer();
        }
        
        function resetZoom() {
            resetView();
        }
        
        // Update zoom display based on viewBox
        function updateZoomDisplayFromViewBox() {
            const svg = document.querySelector('#svgContainer svg');
            if (svg && initialViewBox && initialViewBox.width > 0) {
                const viewBox = svg.viewBox.baseVal;
                if (viewBox.width > 0) {
                    const zoomValue = (initialViewBox.width / viewBox.width).toFixed(1) + 'x';
                    const zoomLevelValueEl = document.getElementById('zoomLevelValue');
                    const zoomValueDisplayEl = document.getElementById('zoomValueDisplay');
                    if (zoomLevelValueEl) zoomLevelValueEl.textContent = zoomValue;
                    if (zoomValueDisplayEl) zoomValueDisplayEl.textContent = zoomValue;
                }
            } else {
                // Fallback to 1.0x if viewBox not initialized
                const zoomLevelValueEl = document.getElementById('zoomLevelValue');
                const zoomValueDisplayEl = document.getElementById('zoomValueDisplay');
                if (zoomLevelValueEl) zoomLevelValueEl.textContent = '1.0x';
                if (zoomValueDisplayEl) zoomValueDisplayEl.textContent = '1.0x';
            }
        }
        
        // Panel dragging functionality
        (function initPanelDrag() {
            let isDragging = false;
            let startX = 0;
            let startY = 0;
            let initialX = 0;
            let initialY = 0;
            
            const panel = document.getElementById('controlPanel');
            const panelHeader = document.getElementById('panelHeader');
            
            if (!panel || !panelHeader) return;
            
            // Get current position
            function getCurrentPosition() {
                const rect = panel.getBoundingClientRect();
                return {
                    x: rect.left,
                    y: rect.top
                };
            }
            
            function dragStart(e) {
                // Don't drag when clicking toggle button or its children
                if (e.target.classList.contains('panel-toggle') || 
                    e.target.closest('.panel-toggle')) {
                    return;
                }
                
                if (e.target === panelHeader || panelHeader.contains(e.target)) {
                    isDragging = true;
                    const pos = getCurrentPosition();
                    initialX = pos.x;
                    initialY = pos.y;
                    startX = e.clientX;
                    startY = e.clientY;
                    e.preventDefault();
                    e.stopPropagation();
                }
            }
            
            function drag(e) {
                if (isDragging) {
                    e.preventDefault();
                    e.stopPropagation();
                    const deltaX = e.clientX - startX;
                    const deltaY = e.clientY - startY;
                    
                    let newX = initialX + deltaX;
                    let newY = initialY + deltaY;
                    
                    // Constrain to viewport
                    const maxX = window.innerWidth - panel.offsetWidth;
                    const maxY = window.innerHeight - panel.offsetHeight;
                    newX = Math.max(0, Math.min(newX, maxX));
                    newY = Math.max(0, Math.min(newY, maxY));
                    
                    panel.style.left = newX + 'px';
                    panel.style.top = newY + 'px';
                    panel.style.right = 'auto';
                }
            }
            
            function dragEnd(e) {
                if (isDragging) {
                    isDragging = false;
                    e.preventDefault();
                    e.stopPropagation();
                }
            }
            
            panelHeader.addEventListener('mousedown', dragStart);
            document.addEventListener('mousemove', drag);
            document.addEventListener('mouseup', dragEnd);
        })();
        
        // Reactivity data embedded in JavaScript
        const reactivityData = {"###);
    
    // Serialize reactivity data as JSON
    html.push_str(&format!("signalType: '{}', positions: {{", signal_type));
    for (idx, data) in reactivity_data.iter().enumerate() {
        if idx > 0 {
            html.push_str(", ");
        }
        html.push_str(&format!(
            "{}: {{base: '{}', reactivity: {:.6}}}",
            data.position, data.base, data.reactivity
        ));
    }
    html.push_str("}};\n");
    
    // METADATA-DRIVEN ARCHITECTURE: Serialize base metadata as JSON
    if let Some(metadata) = base_metadata {
        html.push_str("// Base metadata for metadata-driven visualization\n");
        html.push_str("const RNA_METADATA = [\n");
        for (idx, meta) in metadata.iter().enumerate() {
            if idx > 0 {
                html.push_str(",\n");
            }
            html.push_str(&format!(
                "  {{id: {}, base: '{}', x: {:.2}, y: {:.2}, reactivity: {:.6}}}",
                meta.id, meta.base, meta.x, meta.y, meta.reactivity
            ));
        }
        html.push_str("\n];\n");
    } else {
        html.push_str("// Base metadata not available\n");
        html.push_str("const RNA_METADATA = [];\n");
    }
    
    // Serialize color ranges
    html.push_str("const colorRanges = {\n");
    html.push_str("A: [");
    for (idx, range) in color_ranges.a.iter().enumerate() {
        if idx > 0 { html.push_str(", "); }
        html.push_str(&format!("{{min: {:.2}, max: {:.2}, color: '{}'}}", range.min, range.max, range.color));
    }
    html.push_str("],\n");
    html.push_str("T: [");
    for (idx, range) in color_ranges.t.iter().enumerate() {
        if idx > 0 { html.push_str(", "); }
        html.push_str(&format!("{{min: {:.2}, max: {:.2}, color: '{}'}}", range.min, range.max, range.color));
    }
    html.push_str("],\n");
    html.push_str("C: [");
    for (idx, range) in color_ranges.c.iter().enumerate() {
        if idx > 0 { html.push_str(", "); }
        html.push_str(&format!("{{min: {:.2}, max: {:.2}, color: '{}'}}", range.min, range.max, range.color));
    }
    html.push_str("],\n");
    html.push_str("G: [");
    for (idx, range) in color_ranges.g.iter().enumerate() {
        if idx > 0 { html.push_str(", "); }
        html.push_str(&format!("{{min: {:.2}, max: {:.2}, color: '{}'}}", range.min, range.max, range.color));
    }
    html.push_str("]\n};\n");
    
    // Set initial circle style
    html.push_str(&format!("let circleFilled = {};\n", if initial_circle_filled { "true" } else { "false" }));
    if initial_circle_filled {
        html.push_str("document.getElementById('circleStyle').value = 'filled';\n");
    } else {
        html.push_str("document.getElementById('circleStyle').value = 'hollow';\n");
    }
    
    html.push_str(r###"
        
        // Initialize interactive features
        let currentZoom = 1.0;
        let panStart = { x: 0, y: 0, viewBoxX: 0, viewBoxY: 0 };
        let isPanning = false;
        let highlightedCircles = new Set();
        let initialViewBox = null;
        
        // Tool mode: 'select', 'box', 'hand'
        let currentToolMode = 'select';
        let isSelecting = false;
        let selectionStart = {x: 0, y: 0};
        let selectedBases = new Set();
        
        // Function to initialize SVG interactivity (can be called multiple times)
        function initializeSvgInteractivity() {
            // Get SVG element
            let svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            // Reset selection state
            isSelecting = false;
            selectedBases.clear();
            isPanning = false;
            
            // Hide selection box
            const selectionBox = document.getElementById('selectionBox');
            if (selectionBox) {
                selectionBox.style.display = 'none';
            }
            
            // Remove old event listeners by cloning SVG (removes all listeners)
            const svgParent = svg.parentNode;
            const svgClone = svg.cloneNode(true);
            svgParent.replaceChild(svgClone, svg);
            svg = svgClone;
            
            // Make SVG responsive
            svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
            
            // Store initial viewBox - handle case where viewBox might not be set
            let viewBox = svg.viewBox.baseVal;
            if (viewBox.width === 0 || viewBox.height === 0) {
                // If viewBox is not set, use SVG dimensions
                try {
                    const bbox = svg.getBBox();
                    const svgWidth = parseFloat(svg.getAttribute('width')) || bbox.width;
                    const svgHeight = parseFloat(svg.getAttribute('height')) || bbox.height;
                    if (svgWidth > 0 && svgHeight > 0) {
                        svg.setAttribute('viewBox', `0 0 ${svgWidth} ${svgHeight}`);
                        viewBox = svg.viewBox.baseVal;
                    }
                } catch (e) {
                    console.warn('Could not get SVG dimensions, using defaults');
                    svg.setAttribute('viewBox', '0 0 1000 1000');
                    viewBox = svg.viewBox.baseVal;
                }
            }
            initialViewBox = {
                x: viewBox.x,
                y: viewBox.y,
                width: viewBox.width,
                height: viewBox.height
            };
            
            // Make initialViewBox globally accessible
            window.initialViewBox = initialViewBox;
            
            const svgElement = svg;
            
            // Mouse wheel zoom
            svgElement.addEventListener('wheel', (e) => {
                e.preventDefault();
                const viewBox = svgElement.viewBox.baseVal;
                const delta = e.deltaY > 0 ? 1.1 : 0.9;
                const rect = svgElement.getBoundingClientRect();
                const mouseX = e.clientX - rect.left;
                const mouseY = e.clientY - rect.top;
                
                const x = (mouseX / rect.width) * viewBox.width + viewBox.x;
                const y = (mouseY / rect.height) * viewBox.height + viewBox.y;
                
                viewBox.x = x - (x - viewBox.x) * delta;
                viewBox.y = y - (y - viewBox.y) * delta;
                viewBox.width *= delta;
                viewBox.height *= delta;
                updateZoomDisplayFromViewBox();
            });
            
            // Pan functionality - use spacebar or middle mouse button
            let panMode = false;
            
            // Global key handlers (only add once)
            if (!window.svgKeydownHandler) {
                window.svgKeydownHandler = (e) => {
                    const svg = document.querySelector('#svgContainer svg');
                    if (e.code === 'Space' && e.target.tagName !== 'INPUT' && e.target.tagName !== 'TEXTAREA' && svg) {
                        panMode = true;
                        svg.style.cursor = 'grab';
                        e.preventDefault();
                    }
                };
                window.svgKeyupHandler = (e) => {
                    if (e.code === 'Space') {
                        panMode = false;
                        const svg = document.querySelector('#svgContainer svg');
                        if (svg) svg.style.cursor = '';
                        isPanning = false;
                    }
                };
                document.addEventListener('keydown', window.svgKeydownHandler);
                document.addEventListener('keyup', window.svgKeyupHandler);
            }
            
            // Mouse event handlers - use once flag to prevent duplicates
            const handleMousedown = (e) => {
                if (e.target.closest('.control-panel')) return;
                
                const viewBox = svgElement.viewBox.baseVal;
                
                if (currentToolMode === 'box' && e.button === 0) {
                    // Box select mode
                    isSelecting = true;
                    const rect = svgElement.getBoundingClientRect();
                    selectionStart = {
                        x: e.clientX - rect.left,
                        y: e.clientY - rect.top
                    };
                    const selectionBox = document.getElementById('selectionBox');
                    selectionBox.style.display = 'block';
                    selectionBox.style.left = selectionStart.x + 'px';
                    selectionBox.style.top = selectionStart.y + 'px';
                    selectionBox.style.width = '0px';
                    selectionBox.style.height = '0px';
                    e.preventDefault();
                } else if (currentToolMode === 'hand' && e.button === 0) {
                    // Hand tool for panning
                    isPanning = true;
                    panStart = {
                        x: e.clientX,
                        y: e.clientY
                    };
                    panStart.viewBoxX = viewBox.x;
                    panStart.viewBoxY = viewBox.y;
                    svgElement.style.cursor = 'grabbing';
                    e.preventDefault();
                } else if ((e.button === 1 || (e.button === 0 && panMode)) && !e.target.closest('.control-panel')) {
                    // Middle mouse or spacebar panning
                    isPanning = true;
                    panStart.x = e.clientX;
                    panStart.y = e.clientY;
                    panStart.viewBoxX = viewBox.x;
                    panStart.viewBoxY = viewBox.y;
                    svgElement.style.cursor = 'grabbing';
                    e.preventDefault();
                }
            };
            
            const handleMousemove = (e) => {
                const viewBox = svgElement.viewBox.baseVal;
                
                if (isSelecting) {
                    const currentX = e.clientX;
                    const currentY = e.clientY;
                    const selectionBox = document.getElementById('selectionBox');
                    const width = Math.abs(currentX - selectionStart.x);
                    const height = Math.abs(currentY - selectionStart.y);
                    selectionBox.style.left = Math.min(currentX, selectionStart.x) + 'px';
                    selectionBox.style.top = Math.min(currentY, selectionStart.y) + 'px';
                    selectionBox.style.width = width + 'px';
                    selectionBox.style.height = height + 'px';
                    
                    // Find circles in selection
                    updateSelection();
                } else if (isPanning) {
                    const rect = svgElement.getBoundingClientRect();
                    const scaleX = viewBox.width / rect.width;
                    const scaleY = viewBox.height / rect.height;
                    
                    const deltaX = (panStart.x - e.clientX) * scaleX;
                    const deltaY = (panStart.y - e.clientY) * scaleY;
                    viewBox.x = panStart.viewBoxX + deltaX;
                    viewBox.y = panStart.viewBoxY + deltaY;
                    updateMinimap();
                    updateZoomDisplayFromViewBox();
                    // Reset minimap auto-hide timer on pan
                    resetMinimapHideTimer();
                }
            };
            
            const handleMouseup = (e) => {
                if (isSelecting) {
                    isSelecting = false;
                    document.getElementById('selectionBox').style.display = 'none';
                    if (selectedBases.size > 0) {
                        showBatchEditDialog(e);
                    }
                }
                if (isPanning) {
                    isPanning = false;
                    svgElement.style.cursor = currentToolMode === 'hand' ? 'grab' : '';
                }
            };
            
            // Define event handlers
            const handleWheel = (e) => {
                e.preventDefault();
                const viewBox = svgElement.viewBox.baseVal;
                const delta = e.deltaY > 0 ? 1.1 : 0.9;
                const rect = svgElement.getBoundingClientRect();
                const mouseX = e.clientX - rect.left;
                const mouseY = e.clientY - rect.top;
                
                const x = (mouseX / rect.width) * viewBox.width + viewBox.x;
                const y = (mouseY / rect.height) * viewBox.height + viewBox.y;
                
                viewBox.x = x - (x - viewBox.x) * delta;
                viewBox.y = y - (y - viewBox.y) * delta;
                viewBox.width *= delta;
                viewBox.height *= delta;
                updateZoomDisplayFromViewBox();
                updateMinimap();
                // Reset minimap auto-hide timer on wheel zoom
                resetMinimapHideTimer();
            };
            
            const handleMouseleave = () => {
                if (isPanning) {
                    isPanning = false;
                    if (!panMode) {
                        svgElement.style.cursor = '';
                    }
                }
            };
            
            const handleContextmenu = (e) => {
                e.preventDefault();
            };
            
            // Add event listeners
            svgElement.addEventListener('wheel', handleWheel);
            svgElement.addEventListener('mousedown', handleMousedown);
            svgElement.addEventListener('mousemove', handleMousemove);
            svgElement.addEventListener('mouseup', handleMouseup);
            svgElement.addEventListener('mouseleave', handleMouseleave);
            svgElement.addEventListener('contextmenu', handleContextmenu);
            
            // Add interactivity to circles
            const circles = svgElement.querySelectorAll('circle');
            circles.forEach((circle) => {
                // Extract position from parent's title element or data attribute
                let position = null;
                const parent = circle.parentElement;
                if (parent) {
                    const title = parent.querySelector('title');
                    if (title) {
                        const match = title.textContent.match(/(\d+)/);
                        if (match) {
                            position = parseInt(match[1]);
                        }
                    }
                }
                
                // Also check if data-position already exists
                if (!position && circle.hasAttribute('data-position')) {
                    position = parseInt(circle.getAttribute('data-position'));
                }
                
                if (position && reactivityData.positions && reactivityData.positions[position]) {
                    const data = reactivityData.positions[position];
                    
                    // Add data attributes
                    circle.setAttribute('data-position', position);
                    circle.setAttribute('data-base', data.base);
                    circle.setAttribute('data-reactivity', data.reactivity);
                    
                    // Remove old listeners by cloning
                    const newCircle = circle.cloneNode(true);
                    circle.parentNode.replaceChild(newCircle, circle);
                    
                    // Hover effect
                    newCircle.addEventListener('mouseenter', (e) => {
                        showTooltip(e, position, data);
                        newCircle.style.opacity = '0.8';
                        newCircle.style.cursor = 'pointer';
                    });
                    
                    newCircle.addEventListener('mouseleave', () => {
                        hideTooltip();
                        if (!highlightedCircles.has(position)) {
                            newCircle.style.opacity = '';
                        }
                    });
                    
                    // Click handler - context menu or highlight based on tool mode
                    newCircle.addEventListener('click', (e) => {
                        if (currentToolMode === 'select' && !e.shiftKey && !e.ctrlKey && !e.metaKey) {
                            // Show context menu in select mode
                            e.preventDefault();
                            e.stopPropagation();
                            showContextMenu(e, position, data);
                        } else if (currentToolMode !== 'box' && currentToolMode !== 'hand') {
                            // Toggle highlight in other modes
                            toggleHighlight(newCircle, position);
                        }
                    });
                }
            });
            
            // Initialize selection box
            initSelectionBox();
            initContextMenu();
            
            // Update tool mode cursor
            setToolMode(currentToolMode);
            
            // Ensure batch edit menu exists and is properly initialized
            const batchMenu = document.getElementById('batchEditMenu');
            if (!batchMenu) {
                console.warn('batchEditMenu element not found in DOM');
            }
        }
        
        // Initial call
        const svg = document.querySelector('#svgContainer svg');
        if (svg) {
            initializeSvgInteractivity();
        }
        
        // Context menu functionality
        let currentEditingPosition = null;
        let currentEditingData = null;
        
        function initContextMenu() {
            // Remove old listeners if they exist
            if (window.contextMenuClickHandler) {
                document.removeEventListener('click', window.contextMenuClickHandler);
            }
            
            // Create new handler
            window.contextMenuClickHandler = (e) => {
                const contextMenu = document.getElementById('contextMenu');
                const batchMenu = document.getElementById('batchEditMenu');
                
                if (contextMenu && !e.target.closest('#contextMenu') && !e.target.closest('circle[data-position]')) {
                    contextMenu.style.display = 'none';
                }
                if (batchMenu && !e.target.closest('#batchEditMenu') && !e.target.closest('#selectionBox')) {
                    // Don't close if clicking on selection box
                    if (currentToolMode !== 'box' || !isSelecting) {
                        batchMenu.style.display = 'none';
                    }
                }
            };
            
            document.addEventListener('click', window.contextMenuClickHandler);
        }
        
        function showContextMenu(event, position, data) {
            currentEditingPosition = position;
            currentEditingData = data;
            const menu = document.getElementById('contextMenu');
            menu.style.display = 'block';
            menu.style.left = event.clientX + 'px';
            menu.style.top = event.clientY + 'px';
        }
        
        function editBaseFontSize() {
            if (!currentEditingPosition) return;
            const newSize = prompt('Enter font size (px):', '12');
            if (newSize && !isNaN(parseFloat(newSize)) && parseFloat(newSize) > 0) {
                const text = findTextForPosition(currentEditingPosition);
                if (text) {
                    const sizeValue = parseFloat(newSize) + 'px';
                    text.setAttribute('font-size', sizeValue);
                    text.style.fontSize = sizeValue;
                    saveStateForUndo();
                } else {
                    alert('Text element not found for this position');
                }
            }
            document.getElementById('contextMenu').style.display = 'none';
        }
        
        function editBaseFontColor() {
            if (!currentEditingPosition) return;
            const newColor = prompt('Enter font color (hex, e.g., #000000):', '#000000');
            if (newColor && /^#[0-9A-Fa-f]{6}$/.test(newColor)) {
                const text = findTextForPosition(currentEditingPosition);
                if (text) {
                    text.setAttribute('fill', newColor);
                    text.style.fill = newColor;
                    saveStateForUndo();
                } else {
                    alert('Text element not found for this position');
                }
            } else if (newColor) {
                alert('Invalid color format. Please use hex format (e.g., #000000)');
            }
            document.getElementById('contextMenu').style.display = 'none';
        }
        
        function editBaseCircleRadius() {
            if (!currentEditingPosition) return;
            const newRadius = prompt('Enter circle radius (px):', '2.8');
            if (newRadius && !isNaN(parseFloat(newRadius)) && parseFloat(newRadius) > 0) {
                const circle = document.querySelector(`circle[data-position="${currentEditingPosition}"]`);
                if (circle) {
                    const radiusValue = parseFloat(newRadius);
                    circle.setAttribute('r', radiusValue);
                    // Also update style if radius is set via style
                    if (circle.style.r) {
                        circle.style.r = radiusValue + 'px';
                    }
                    saveStateForUndo();
                } else {
                    alert('Circle element not found for this position');
                }
            } else if (newRadius) {
                alert('Invalid radius. Please enter a positive number');
            }
            document.getElementById('contextMenu').style.display = 'none';
        }
        
        function editBaseCircleColor() {
            if (!currentEditingPosition) return;
            const newColor = prompt('Enter circle color (hex, e.g., #ff0000):', '#000000');
            if (newColor && /^#[0-9A-Fa-f]{6}$/.test(newColor)) {
                const circle = document.querySelector(`circle[data-position="${currentEditingPosition}"]`);
                if (circle) {
                    // Get current style to preserve stroke-width and other properties
                    const currentFill = circle.style.fill || window.getComputedStyle(circle).fill;
                    const currentStroke = circle.style.stroke || window.getComputedStyle(circle).stroke;
                    const currentStrokeWidth = circle.style.strokeWidth || window.getComputedStyle(circle).strokeWidth;
                    
                    // Determine if circle is filled or hollow based on current fill
                    const isFilled = currentFill && currentFill !== 'none' && currentFill !== 'transparent';
                    
                    // Update fill and stroke using DOM style object for reliable updates
                    if (isFilled) {
                        // Filled circle: set both fill and stroke to new color
                        circle.style.fill = newColor;
                        circle.style.stroke = newColor;
                    } else {
                        // Hollow circle: keep fill as none, update stroke
                        circle.style.fill = 'none';
                        circle.style.stroke = newColor;
                    }
                    
                    // Preserve stroke-width if it exists
                    if (currentStrokeWidth && currentStrokeWidth !== '0px') {
                        circle.style.strokeWidth = currentStrokeWidth;
                    }
                    
                    // Also update attributes for consistency
                    circle.setAttribute('fill', isFilled ? newColor : 'none');
                    circle.setAttribute('stroke', newColor);
                    
                    saveStateForUndo();
                } else {
                    alert('Circle element not found for this position');
                }
            } else if (newColor) {
                alert('Invalid color format. Please use hex format (e.g., #ff0000)');
            }
            document.getElementById('contextMenu').style.display = 'none';
        }
        
        function resetBaseStyle() {
            if (!currentEditingPosition) return;
            const text = findTextForPosition(currentEditingPosition);
            const circle = document.querySelector(`circle[data-position="${currentEditingPosition}"]`);
            
            if (text) {
                text.setAttribute('font-size', '12px');
                text.setAttribute('fill', '#000000');
                text.style.fontSize = '12px';
                text.style.fill = '#000000';
            }
            
            if (circle) {
                // Reset radius
                circle.setAttribute('r', '2.8');
                if (circle.style.r) {
                    circle.style.r = '2.8px';
                }
                
                // Note: We don't reset circle color here as the original color
                // depends on reactivity data. User should use batch reset if needed.
            }
            
            saveStateForUndo();
            document.getElementById('contextMenu').style.display = 'none';
        }
        
        function findTextForPosition(position) {
            const positionStr = String(position);
            
            // Method 1: Try to find text element in the same group as circle
            const circles = document.querySelectorAll(`circle[data-position="${position}"]`);
            if (circles.length > 0) {
                const circle = circles[0];
                const parent = circle.parentElement;
                if (parent) {
                    const textInGroup = parent.querySelector('text');
                    if (textInGroup) {
                        return textInGroup;
                    }
                }
            }
            
            // Method 2: Find text element by title tag matching the position
            // Text elements from template have <title>position (position.label in template: ...)</title>
            // or <title>position</title> before them
            const allGroups = document.querySelectorAll('#svgContainer svg g');
            for (let group of allGroups) {
                const title = group.querySelector('title');
                if (title) {
                    const titleText = title.textContent.trim();
                    // Match title that starts with the position number
                    // Examples: "790", "790 (position.label in template: 788.U)"
                    if (titleText === positionStr || titleText.startsWith(positionStr + ' ')) {
                        const text = group.querySelector('text');
                        if (text) {
                            return text;
                        }
                    }
                }
            }
            
            // Method 3: Find text element by searching nearby groups (text might be in adjacent group)
            if (circles.length > 0) {
                const circle = circles[0];
                const circleGroup = circle.parentElement;
                if (circleGroup) {
                    // Check previous and next sibling groups
                    let sibling = circleGroup.previousElementSibling;
                    if (sibling) {
                        const text = sibling.querySelector('text');
                        if (text) return text;
                    }
                    sibling = circleGroup.nextElementSibling;
                    if (sibling) {
                        const text = sibling.querySelector('text');
                        if (text) return text;
                    }
                }
            }
            
            return null;
        }
        
        function setToolMode(mode) {
            currentToolMode = mode;
            // Update button states
            document.querySelectorAll('.tool-btn').forEach(btn => {
                btn.classList.remove('active');
            });
            const activeBtn = document.querySelector(`[data-tool="${mode}"]`);
            if (activeBtn) {
                activeBtn.classList.add('active');
            }
            // Update cursor
            const svg = document.querySelector('#svgContainer svg');
            if (svg) {
                if (mode === 'hand') {
                    svg.style.cursor = 'grab';
                } else if (mode === 'box') {
                    svg.style.cursor = 'crosshair';
                } else {
                    svg.style.cursor = 'default';
                }
            }
        }
        
        function initSelectionBox() {
            // Selection box initialization is now handled in initializeSvgInteractivity
            // This function is kept for compatibility but does nothing
        }
        
        function updateSelection() {
            const selectionBox = document.getElementById('selectionBox');
            const boxRect = selectionBox.getBoundingClientRect();
            
            selectedBases.clear();
            const circles = document.querySelectorAll('#svgContainer svg circle[data-position]');
            circles.forEach(circle => {
                const circleRect = circle.getBoundingClientRect();
                if (circleRect.left >= boxRect.left && circleRect.right <= boxRect.right &&
                    circleRect.top >= boxRect.top && circleRect.bottom <= boxRect.bottom) {
                    const position = circle.getAttribute('data-position');
                    selectedBases.add(position);
                    circle.classList.add('highlighted');
                } else {
                    if (!highlightedCircles.has(parseInt(circle.getAttribute('data-position')))) {
                        circle.classList.remove('highlighted');
                    }
                }
            });
        }
        
        function showBatchEditDialog(event) {
            if (selectedBases.size === 0) return;
            
            const menu = document.getElementById('batchEditMenu');
            if (!menu) {
                console.error('batchEditMenu element not found');
                return;
            }
            
            const countEl = document.getElementById('batchEditCount');
            if (countEl) {
                countEl.textContent = selectedBases.size;
            }
            
            // Position menu at selection center or event position
            if (event) {
                menu.style.left = event.clientX + 'px';
                menu.style.top = event.clientY + 'px';
            } else {
                // Calculate center of selection
                let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
                selectedBases.forEach(position => {
                    const circle = document.querySelector(`circle[data-position="${position}"]`);
                    if (circle) {
                        const rect = circle.getBoundingClientRect();
                        minX = Math.min(minX, rect.left + rect.width / 2);
                        minY = Math.min(minY, rect.top + rect.height / 2);
                        maxX = Math.max(maxX, rect.left + rect.width / 2);
                        maxY = Math.max(maxY, rect.top + rect.height / 2);
                    }
                });
                if (minX !== Infinity) {
                    menu.style.left = ((minX + maxX) / 2) + 'px';
                    menu.style.top = ((minY + maxY) / 2) + 'px';
                }
            }
            
            menu.style.display = 'block';
            menu.style.zIndex = '10001';
            
            // Ensure menu is visible and properly positioned
            setTimeout(() => {
                const rect = menu.getBoundingClientRect();
                if (rect.right > window.innerWidth) {
                    menu.style.left = (window.innerWidth - rect.width - 10) + 'px';
                }
                if (rect.bottom > window.innerHeight) {
                    menu.style.top = (window.innerHeight - rect.height - 10) + 'px';
                }
            }, 0);
        }
        
        function batchEditFontSize() {
            if (selectedBases.size === 0) return;
            const newSize = prompt('Enter font size (px):', '12');
            if (newSize) {
                selectedBases.forEach(position => {
                    const text = findTextForPosition(position);
                    if (text) {
                        text.setAttribute('font-size', newSize + 'px');
                        text.style.fontSize = newSize + 'px';
                    }
                });
                saveStateForUndo();
            }
            closeBatchEditMenu();
        }
        
        function batchEditFontColor() {
            if (selectedBases.size === 0) return;
            const newColor = prompt('Enter font color (hex):', '#000000');
            if (newColor) {
                selectedBases.forEach(position => {
                    const text = findTextForPosition(position);
                    if (text) {
                        text.setAttribute('fill', newColor);
                        text.style.fill = newColor;
                    }
                });
                saveStateForUndo();
            }
            closeBatchEditMenu();
        }
        
        function batchEditCircleRadius() {
            if (selectedBases.size === 0) return;
            const newRadius = prompt('Enter circle radius (px):', '2.8');
            if (newRadius) {
                selectedBases.forEach(position => {
                    const circle = document.querySelector(`circle[data-position="${position}"]`);
                    if (circle) {
                        circle.setAttribute('r', newRadius);
                    }
                });
                saveStateForUndo();
            }
            closeBatchEditMenu();
        }
        
        function batchEditCircleColor() {
            if (selectedBases.size === 0) return;
            const newColor = prompt('Enter circle color (hex):', '#000000');
            if (newColor) {
                selectedBases.forEach(position => {
                    const circle = document.querySelector(`circle[data-position="${position}"]`);
                    if (circle) {
                        const currentStyle = circle.getAttribute('style') || '';
                        const newStyle = currentStyle.replace(/fill:[^;]+;?/g, '').replace(/stroke:[^;]+;?/g, '');
                        circle.setAttribute('style', newStyle + `fill: ${newColor}; stroke: ${newColor};`);
                    }
                });
                saveStateForUndo();
            }
            closeBatchEditMenu();
        }
        
        function batchResetStyle() {
            if (selectedBases.size === 0) return;
            selectedBases.forEach(position => {
                const text = findTextForPosition(position);
                const circle = document.querySelector(`circle[data-position="${position}"]`);
                if (text) {
                    text.setAttribute('font-size', '12px');
                    text.setAttribute('fill', '#000000');
                }
                if (circle) {
                    circle.setAttribute('r', '2.8');
                }
            });
            saveStateForUndo();
            closeBatchEditMenu();
        }
        
        function closeBatchEditMenu() {
            const menu = document.getElementById('batchEditMenu');
            menu.style.display = 'none';
            selectedBases.forEach(position => {
                const circle = document.querySelector(`circle[data-position="${position}"]`);
                if (circle && !highlightedCircles.has(parseInt(position))) {
                    circle.classList.remove('highlighted');
                }
            });
            selectedBases.clear();
        }
        
        // Tooltip functions
        function showTooltip(event, position, data) {
            const tooltip = document.getElementById('tooltip');
            tooltip.innerHTML = `
                <strong>Position:</strong> ${position}<br>
                <strong>Base:</strong> ${data.base}<br>
                <strong>Reactivity:</strong> ${data.reactivity.toFixed(4)}<br>
                <strong>Signal:</strong> ${reactivityData.signalType}
            `;
            tooltip.classList.add('visible');
            updateTooltipPosition(event);
        }
        
        function hideTooltip() {
            const tooltip = document.getElementById('tooltip');
            tooltip.classList.remove('visible');
        }
        
        function updateTooltipPosition(event) {
            const tooltip = document.getElementById('tooltip');
            tooltip.style.left = (event.clientX + 10) + 'px';
            tooltip.style.top = (event.clientY + 10) + 'px';
        }
        
        // Highlight functions
        function toggleHighlight(circle, position) {
            if (highlightedCircles.has(position)) {
                circle.classList.remove('highlighted');
                highlightedCircles.delete(position);
            } else {
                circle.classList.add('highlighted');
                highlightedCircles.add(position);
            }
        }
        
        // Determine color based on reactivity and color ranges
        function determineColor(reactivity, base) {
            // CRITICAL FIX: Normalize T to U because colorRanges uses U as the key
            if (base) {
                base = base.toUpperCase().trim();
                // Normalize T to U (colorRanges uses U as the key, not T)
                if (base === 'T') {
                    base = 'U';
                }
            }
            
            // Get color ranges for this specific base (now normalized to U if it was T)
            const ranges = colorRanges[base];
            if (!ranges) {
                // Fallback: try T if U not found (for legacy compatibility)
                if (base === 'U' && colorRanges['T']) {
                    return determineColorFromRanges(reactivity, colorRanges['T']);
                }
                return 'rgba(169, 169, 169, 0.8)';
            }
            
            return determineColorFromRanges(reactivity, ranges);
        }
        
        // Helper function to determine color from ranges
        function determineColorFromRanges(reactivity, ranges) {
            
            // Handle invalid reactivity values
            if (reactivity === null || reactivity === undefined || isNaN(reactivity) || reactivity < 0) {
                return 'rgba(169, 169, 169, 0.8)';
            }
            
            // Check ranges in order, using proper boundary conditions
            // For ranges: [min, max) except the last one which includes max
            for (let i = 0; i < ranges.length; i++) {
                const range = ranges[i];
                const isLastRange = (i === ranges.length - 1);
                
                // For all ranges except the last: min <= reactivity < max
                // For the last range: min <= reactivity <= max
                if (isLastRange) {
                    if (reactivity >= range.min && reactivity <= range.max) {
                        return range.color;
                    }
                } else {
                    if (reactivity >= range.min && reactivity < range.max) {
                        return range.color;
                    }
                }
            }
            
            // If reactivity is outside all ranges, return gray
            return 'rgba(169, 169, 169, 0.8)';
        }
        
        // Convert hex color to rgba
        function hexToRgba(hex, alpha = 0.8) {
            const r = parseInt(hex.slice(1, 3), 16);
            const g = parseInt(hex.slice(3, 5), 16);
            const b = parseInt(hex.slice(5, 7), 16);
            return `rgba(${r}, ${g}, ${b}, ${alpha})`;
        }
        
        // Helper function to get color value from Spectrum or regular input
        function getColorValue(inputId) {
            const input = document.getElementById(inputId);
            if (!input) return '#000000';
            
            // First try to get value from input directly (most reliable)
            if (input.value) {
                // If it's already a valid hex color, return it
                if (/^#[0-9A-Fa-f]{6}$/.test(input.value)) {
                    return input.value;
                }
            }
            
            // If Spectrum is initialized, try to get the color value
            if (typeof $ !== 'undefined') {
                try {
                    // Check if Spectrum is actually initialized on this element
                    const spectrumInstance = $(input).data('spectrum');
                    if (!spectrumInstance) {
                        // Spectrum not initialized yet, use input value
                        return input.value || '#000000';
                    }
                    
                    const color = $(input).spectrum('get');
                    if (!color) {
                        return input.value || '#000000';
                    }
                    
                    // Check if color is a Spectrum color object
                    if (color && typeof color === 'object') {
                        // First check if it has toHexString method
                        if (typeof color.toHexString === 'function') {
                            try {
                                return color.toHexString();
                            } catch (e) {
                                // If toHexString fails, try toString
                            }
                        }
                        // Some Spectrum versions return color as object with toString
                        if (typeof color.toString === 'function') {
                            try {
                                const colorStr = color.toString();
                                if (colorStr && (colorStr.startsWith('#') || /^rgb|rgba|hsl|hsla/i.test(colorStr))) {
                                    return colorStr;
                                }
                            } catch (e) {
                                // If toString fails, fall back to input value
                            }
                        }
                        // If object has properties that look like color values
                        if (color.r !== undefined && color.g !== undefined && color.b !== undefined) {
                            const r = Math.round(color.r);
                            const g = Math.round(color.g);
                            const b = Math.round(color.b);
                            return `#${r.toString(16).padStart(2, '0')}${g.toString(16).padStart(2, '0')}${b.toString(16).padStart(2, '0')}`;
                        }
                    }
                    // If color is already a string, return it
                    if (typeof color === 'string') {
                        if (color.startsWith('#') || /^rgb|rgba|hsl|hsla/i.test(color)) {
                            return color;
                        }
                    }
                    // If we get here and color is an object but none of the methods worked, use input value
                    return input.value || '#000000';
                } catch (e) {
                    // Silently fail and use input value
                    console.warn('Error getting Spectrum color for', inputId, ':', e);
                    return input.value || '#000000';
                }
            }
            
            // Fallback to input value or default
            return input.value || '#000000';
        }
        
        // Toggle unified mode
        function toggleUnifiedMode() {
            const unified = document.getElementById('unifiedColorRanges').checked;
            const unifiedControls = document.getElementById('unified-controls');
            const individualControls = document.getElementById('individual-controls');
            
            if (unified) {
                unifiedControls.style.display = 'block';
                individualControls.style.display = 'none';
                // Sync unified controls with A's values
                const range1 = parseFloat(document.getElementById('rangeA1').value);
                const range2 = parseFloat(document.getElementById('rangeA2').value);
                document.getElementById('unifiedRange1').value = range1;
                document.getElementById('unifiedRange2').value = range2;
                document.getElementById('unifiedRange1Value').textContent = range1.toFixed(2);
                document.getElementById('unifiedRange2Value').textContent = range2.toFixed(2);
                document.getElementById('unifiedColor1').value = getColorValue('colorA1');
                document.getElementById('unifiedColor2').value = getColorValue('colorA2');
                document.getElementById('unifiedColor3').value = getColorValue('colorA3');
                if (typeof $ !== 'undefined') {
                    $('#unifiedColor1').spectrum('set', getColorValue('colorA1'));
                    $('#unifiedColor2').spectrum('set', getColorValue('colorA2'));
                    $('#unifiedColor3').spectrum('set', getColorValue('colorA3'));
                }
                // CRITICAL: Only update color ranges when switching TO unified mode
                // This ensures all bases get the unified settings
                updateColorRanges();
                // Update all bases using verified metadata-driven approach
                if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                    ['A', 'U', 'C', 'G'].forEach(base => {
                        RNAView.updateByBase(base);
                    });
                } else {
                    applyColorAndVisibility();
                }
            } else {
                unifiedControls.style.display = 'none';
                individualControls.style.display = 'block';
                // CRITICAL: When switching to individual mode, update color ranges from individual controls
                // This ensures each base uses its own settings
            updateColorRanges();
                // Update all bases using verified metadata-driven approach
                if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                    ['A', 'U', 'C', 'G'].forEach(base => {
                        RNAView.updateByBase(base);
                    });
                } else {
                    applyColorAndVisibility();
                }
            }
            updateDistributionPreviews();
            updateLegendValues();
        }
        
        // Draw reactivity distribution preview
        function drawDistributionPreview(base, canvasId, threshold1, threshold2, color1, color2, color3) {
            const canvas = document.getElementById(canvasId);
            if (!canvas) {
                console.warn('Canvas not found:', canvasId);
                return;
            }
            
            console.log('Drawing distribution for base:', base, 'canvasId:', canvasId, 'thresholds:', threshold1, threshold2);
            
            // Set canvas size to match container
            const container = canvas.parentElement;
            let width = 300;
            let height = 60;
            
            if (container) {
                const rect = container.getBoundingClientRect();
                width = rect.width || 300;
                height = rect.height || 60;
                // Use actual pixel dimensions for crisp rendering
                const dpr = window.devicePixelRatio || 1;
                canvas.width = width * dpr;
                canvas.height = height * dpr;
            } else {
                canvas.width = 300;
                canvas.height = 60;
            }
            
            const ctx = canvas.getContext('2d', { alpha: true });
            if (container && window.devicePixelRatio) {
                ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
            }
            
            // Clear canvas and fill with white background
            ctx.clearRect(0, 0, width, height);
            ctx.fillStyle = 'white';
            ctx.fillRect(0, 0, width, height);
            
            // Get reactivity values for this base
            const reactivities = [];
            if (reactivityData && reactivityData.positions) {
                for (const [pos, data] of Object.entries(reactivityData.positions)) {
                    if (data.base === base) {
                        reactivities.push(data.reactivity);
                    }
                }
            }
            
            console.log('Found', reactivities.length, 'reactivity values for base', base);
            
            if (reactivities.length === 0) {
                // No data, draw white background with text
                ctx.fillStyle = 'white';
                ctx.fillRect(0, 0, width, height);
                ctx.fillStyle = '#999';
                ctx.font = '10px Arial';
                ctx.textAlign = 'center';
                ctx.fillText('No data', width / 2, height / 2);
                return;
            }
            
            // Find data range for normalization
            const minVal = Math.min(...reactivities);
            const maxVal = Math.max(...reactivities);
            const dataRange = maxVal - minVal;
            
            // Create histogram bins (50 bins)
            const bins = 50;
            const histogram = new Array(bins).fill(0);
            
            reactivities.forEach(r => {
                // Normalize to 0-1 range if needed
                let normalizedR = r;
                if (dataRange > 0) {
                    normalizedR = (r - minVal) / dataRange;
                }
                // Clamp to 0-1
                normalizedR = Math.max(0, Math.min(1, normalizedR));
                
                const binIndex = Math.min(Math.floor(normalizedR * bins), bins - 1);
                histogram[binIndex]++;
            });
            
            const maxCount = Math.max(...histogram, 1);
            const binWidth = width / bins;
            
            // Draw histogram bars with colors based on thresholds
            histogram.forEach((count, index) => {
                const x = index * binWidth;
                const barHeight = (count / maxCount) * height;
                
                // Determine color based on bin position (normalized to 0-1)
                const binValue = (index + 0.5) / bins;
                let color;
                if (binValue < threshold1) {
                    color = color1;
                } else if (binValue < threshold2) {
                    color = color2;
                } else {
                    color = color3;
                }
                
                ctx.fillStyle = color;
                ctx.fillRect(x, height - barHeight, binWidth, barHeight);
            });
            
            // Draw threshold lines
            ctx.strokeStyle = '#333';
            ctx.lineWidth = 1;
            const threshold1X = threshold1 * width;
            const threshold2X = threshold2 * width;
            ctx.beginPath();
            ctx.moveTo(threshold1X, 0);
            ctx.lineTo(threshold1X, height);
            ctx.stroke();
            ctx.beginPath();
            ctx.moveTo(threshold2X, 0);
            ctx.lineTo(threshold2X, height);
            ctx.stroke();
        }
        
        // Debounce timer for distribution preview updates
        let distributionPreviewTimeout = null;
        
        // Update distribution previews with debouncing and optional base filter
        function updateDistributionPreviews(baseFilter = null) {
            // Clear existing timeout
            if (distributionPreviewTimeout) {
                clearTimeout(distributionPreviewTimeout);
            }
            
            // Debounce: only update after 150ms of no changes
            distributionPreviewTimeout = setTimeout(() => {
            const unifiedEl = document.getElementById('unifiedColorRanges');
            if (!unifiedEl) {
                console.warn('unifiedColorRanges element not found');
                return;
            }
            const unified = unifiedEl.checked;
            
            if (unified) {
                // Update unified distribution
                const range1 = parseFloat(document.getElementById('unifiedRange1').value);
                const range2 = parseFloat(document.getElementById('unifiedRange2').value);
                const color1 = getColorValue('unifiedColor1');
                const color2 = getColorValue('unifiedColor2');
                const color3 = getColorValue('unifiedColor3');
                
                // Draw for all bases combined
                const allReactivities = [];
                if (reactivityData && reactivityData.positions) {
                    for (const [pos, data] of Object.entries(reactivityData.positions)) {
                        allReactivities.push(data.reactivity);
                    }
                }
                
                const canvas = document.getElementById('unifiedDistributionCanvas');
                if (canvas) {
                    // Set canvas size to match container
                    const container = canvas.parentElement;
                    let width, height;
                    if (container) {
                        const rect = container.getBoundingClientRect();
                        width = rect.width || 300;
                        height = rect.height || 60;
                        // Use actual pixel dimensions
                        const dpr = window.devicePixelRatio || 1;
                        canvas.width = width * dpr;
                        canvas.height = height * dpr;
                    } else {
                        width = 300;
                        height = 60;
                        canvas.width = 300;
                        canvas.height = 60;
                    }
                    
                    const ctx = canvas.getContext('2d', { alpha: true });
                    if (container && window.devicePixelRatio) {
                        ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
                    }
                    // Clear canvas and fill with white background
                    ctx.clearRect(0, 0, width, height);
                    ctx.fillStyle = 'white';
                    ctx.fillRect(0, 0, width, height);
                    
                    if (allReactivities.length > 0) {
                        // Find data range for normalization
                        const minVal = Math.min(...allReactivities);
                        const maxVal = Math.max(...allReactivities);
                        const dataRange = maxVal - minVal;
                        
                        const bins = 50;
                        const histogram = new Array(bins).fill(0);
                        
                        allReactivities.forEach(r => {
                            // Normalize to 0-1 range if needed
                            let normalizedR = r;
                            if (dataRange > 0) {
                                normalizedR = (r - minVal) / dataRange;
                            }
                            // Clamp to 0-1
                            normalizedR = Math.max(0, Math.min(1, normalizedR));
                            
                            const binIndex = Math.min(Math.floor(normalizedR * bins), bins - 1);
                            histogram[binIndex]++;
                        });
                        
                        const maxCount = Math.max(...histogram, 1);
                        const binWidth = width / bins;
                        
                        histogram.forEach((count, index) => {
                            const x = index * binWidth;
                            const barHeight = (count / maxCount) * height;
                            const binValue = (index + 0.5) / bins;
                            
                            let color;
                            if (binValue < range1) {
                                color = color1;
                            } else if (binValue < range2) {
                                color = color2;
                            } else {
                                color = color3;
                            }
                            
                            ctx.fillStyle = color;
                            ctx.fillRect(x, height - barHeight, binWidth, barHeight);
                        });
                        
                        // Draw threshold lines
                        ctx.strokeStyle = '#333';
                        ctx.lineWidth = 1;
                        const threshold1X = range1 * width;
                        const threshold2X = range2 * width;
                        ctx.beginPath();
                        ctx.moveTo(threshold1X, 0);
                        ctx.lineTo(threshold1X, height);
                        ctx.stroke();
                        ctx.beginPath();
                        ctx.moveTo(threshold2X, 0);
                        ctx.lineTo(threshold2X, height);
                        ctx.stroke();
                    } else {
                        // No data, draw white background with text
                        ctx.fillStyle = 'white';
                        ctx.fillRect(0, 0, width, height);
                        ctx.fillStyle = '#999';
                        ctx.font = '10px Arial';
                        ctx.textAlign = 'center';
                        ctx.fillText('No data', width / 2, height / 2);
                    }
                }
            } else {
                // Update individual base distributions
                    // OPTIMIZATION: Only update specified base if filter provided, otherwise update all
                    const basesToUpdate = baseFilter ? [baseFilter] : ['A', 'T', 'C', 'G'];
                    
                    basesToUpdate.forEach(base => {
                    const range1El = document.getElementById(`range${base}1`);
                    const range2El = document.getElementById(`range${base}2`);
                    if (!range1El || !range2El) {
                        console.warn(`Missing range elements for base ${base}`);
                        return;
                    }
                    const range1 = parseFloat(range1El.value);
                    const range2 = parseFloat(range2El.value);
                    
                    // Get colors with error handling
                    let color1, color2, color3;
                    try {
                        color1 = getColorValue(`color${base}1`);
                        color2 = getColorValue(`color${base}2`);
                        color3 = getColorValue(`color${base}3`);
                    } catch (e) {
                        console.warn(`Error getting colors for base ${base}:`, e);
                        // Use default colors
                        color1 = '#f5f5f5';
                        color2 = '#ffa500';
                        color3 = '#ff0000';
                    }
                    
                    drawDistributionPreview(base, `distributionCanvas${base}`, range1, range2, color1, color2, color3);
                });
            }
            }, 150); // 150ms debounce delay
        }
        
        // Update color ranges from controls
        // Update color ranges for a specific base only (for atomic updates)
        function updateColorRangeForBase(base) {
            // Normalize: T -> U
            const normalizedBase = base === 'T' ? 'U' : base;
            
            // Helper function to get and validate range values
            const getRange = (base, idx) => {
                const value = parseFloat(document.getElementById(`range${base}${idx}`).value);
                return Math.max(0.0, Math.min(value, 1.0));
            };
            
            const range1 = getRange(base, 1);
            const range2 = Math.max(range1, getRange(base, 2));
            
            // CRITICAL: Only update the color range for this specific base
            const newRanges = [
                {min: 0.0, max: range1, color: hexToRgba(getColorValue(`color${base}1`))},
                {min: range1, max: range2, color: hexToRgba(getColorValue(`color${base}2`))},
                {min: range2, max: 1.0, color: hexToRgba(getColorValue(`color${base}3`))}
            ];
            
            // CRITICAL FIX: Use deep copy instead of reference assignment to prevent cross-contamination
            colorRanges[normalizedBase] = newRanges.map(range => ({
                min: range.min,
                max: range.max,
                color: range.color
            }));
            
            // Also set for T if this is T/U (for legacy compatibility) - use deep copy
            if (base === 'T') {
                colorRanges['T'] = newRanges.map(range => ({
                    min: range.min,
                    max: range.max,
                    color: range.color
                }));
            }
        }
        
        function updateColorRanges() {
            const unified = document.getElementById('unifiedColorRanges').checked;
            
            if (unified) {
                // Use unified ranges and colors for all bases
                const range1 = parseFloat(document.getElementById('unifiedRange1').value);
                const range2 = parseFloat(document.getElementById('unifiedRange2').value);
                const color1 = hexToRgba(getColorValue('unifiedColor1'));
                const color2 = hexToRgba(getColorValue('unifiedColor2'));
                const color3 = hexToRgba(getColorValue('unifiedColor3'));
                
                // Ensure ranges are valid and ordered
                const validRange1 = Math.max(0.0, Math.min(range1, 1.0));
                const validRange2 = Math.max(validRange1, Math.min(range2, 1.0));
                
                // Update color ranges for all bases (normalize T to U)
                ['A', 'U', 'C', 'G'].forEach(base => {
                    colorRanges[base] = [
                        {min: 0.0, max: validRange1, color: color1},
                        {min: validRange1, max: validRange2, color: color2},
                        {min: validRange2, max: 1.0, color: color3}
                    ];
                });
                // Also update T for legacy compatibility - use deep copy
                colorRanges['T'] = colorRanges['U'].map(range => ({
                    min: range.min,
                    max: range.max,
                    color: range.color
                }));
            } else {
                // Use individual ranges and colors - update ALL bases
                // Helper function to get and validate range values
                const getRange = (base, idx) => {
                    const value = parseFloat(document.getElementById(`range${base}${idx}`).value);
                    return Math.max(0.0, Math.min(value, 1.0));
                };
                
                ['A', 'T', 'C', 'G'].forEach(base => {
                    const range1 = getRange(base, 1);
                    const range2 = Math.max(range1, getRange(base, 2));
                    
                    const normalizedBase = base === 'T' ? 'U' : base;
                    colorRanges[normalizedBase] = [
                        {min: 0.0, max: range1, color: hexToRgba(getColorValue(`color${base}1`))},
                        {min: range1, max: range2, color: hexToRgba(getColorValue(`color${base}2`))},
                        {min: range2, max: 1.0, color: hexToRgba(getColorValue(`color${base}3`))}
                    ];
                    // Also set for original base for legacy compatibility - use deep copy
                    if (base === 'T') {
                        colorRanges['T'] = colorRanges['U'].map(range => ({
                            min: range.min,
                            max: range.max,
                            color: range.color
                        }));
                    }
                });
            }
            
            updateDistributionPreviews();
            updateLegendValues();
            // Use verified metadata-driven approach
            applyColorAndVisibility();
        }
        
        // Update legend values based on current color ranges
        // baseFilter: optional parameter to update only a specific base (e.g., 'A', 'T', 'C', 'G')
        function updateLegendValues(baseFilter = null) {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            const unified = document.getElementById('unifiedColorRanges').checked;
            
            // CRITICAL FIX: Only update specified base if filter provided, otherwise update all
            const basesToUpdate = baseFilter ? [baseFilter] : ['A', 'T', 'C', 'G'];
            
            basesToUpdate.forEach(base => {
                // Find legend group for this base
                const legendGroup = svg.querySelector(`#legend-${base}, g[data-base="${base}"].legend-group`);
                if (!legendGroup) return;
                
                // CRITICAL FIX: Handle T/U normalization - T should use U's color ranges
                const normalizedBase = base === 'T' ? 'U' : base;
                const ranges = colorRanges[normalizedBase] || colorRanges[base];
                if (!ranges || ranges.length < 3) return;
                
                // Find legend items (excluding the -1 item)
                const itemsGroup = legendGroup.querySelector('.legend-items');
                if (!itemsGroup) return;
                
                const items = Array.from(itemsGroup.querySelectorAll('.legend-item[data-min]'));
                
                // Update each range item (skip the -1 item which has data-value instead of data-min)
                items.forEach((item, idx) => {
                    if (idx < ranges.length) {
                        const range = ranges[idx];
                        const text = item.querySelector('text.legend-text, text');
                        
                        if (text && range) {
                            // Update text content
                            text.textContent = `${range.min.toFixed(2)} - ${range.max.toFixed(2)}`;
                            
                            // Update data attributes
                            item.setAttribute('data-min', range.min.toFixed(2));
                            item.setAttribute('data-max', range.max.toFixed(2));
                        }
                    }
                });
            });
        }
        
        // ============================================================================
        // COMPREHENSIVE VERIFICATION ARCHITECTURE
        // ============================================================================
        // Layer 1: SVG Template Reference (Read-only verification source)
        // Layer 2: Metadata Structure (Single source of truth)
        // Layer 3: Independent Visualization Layer (Metadata-driven rendering)
        // Layer 4: Control Panel (Validated updates)
        // ============================================================================
        
        // Layer 1: SVG Template Validator - Extract reference data from original SVG
        const SVGTemplateValidator = {
            // Extract position and base information from original SVG template (read-only)
            extractTemplateReference() {
                const svg = document.querySelector('#svgContainer svg');
                if (!svg) return null;
                
                const reference = {
                    positions: new Map(),  // position -> {base, x, y}
                    circles: new Map()      // node_id -> {position, base, x, y}
                };
                
                // Extract from all circles with node IDs
                const allCircles = svg.querySelectorAll('circle[id^="node_"]');
                allCircles.forEach(circle => {
                    const id = circle.getAttribute('id');
                    const position = parseInt(circle.getAttribute('data-position'));
                    const base = circle.getAttribute('data-base');
                    const x = parseFloat(circle.getAttribute('cx'));
                    const y = parseFloat(circle.getAttribute('cy'));
                    
                    if (id && !isNaN(position) && base) {
                        const normalizedBase = base.toUpperCase() === 'T' ? 'U' : base.toUpperCase();
                        reference.positions.set(position, {base: normalizedBase, x, y});
                        reference.circles.set(id, {position, base: normalizedBase, x, y});
                    }
                });
                
                console.log('SVG Template Reference extracted:', {
                    totalPositions: reference.positions.size,
                    totalCircles: reference.circles.size
                });
                
                return reference;
            },
            
            // Verify metadata against SVG template
            verifyMetadata(metadata, templateRef) {
                if (!templateRef) {
                    console.warn('No template reference available for verification');
                    return {valid: false, errors: ['No template reference']};
                }
                
                const errors = [];
                const warnings = [];
                let verifiedCount = 0;
                
                metadata.forEach(meta => {
                    const nodeId = `node_${meta.id}`;
                    const templateCircle = templateRef.circles.get(nodeId);
                    const templatePos = templateRef.positions.get(meta.id);
                    
                    // Verify circle exists in template
                    if (!templateCircle) {
                        errors.push(`Metadata ID ${meta.id} (node_${meta.id}) not found in SVG template`);
                        return;
                    }
                    
                    // Verify position matches
                    if (templateCircle.position !== meta.id) {
                        errors.push(`Position mismatch for ID ${meta.id}: template=${templateCircle.position}, metadata=${meta.id}`);
                        return;
                    }
                    
                    // Verify base matches (with T/U normalization)
                    const metaBase = meta.base.toUpperCase() === 'T' ? 'U' : meta.base.toUpperCase();
                    if (templateCircle.base !== metaBase) {
                        errors.push(`Base mismatch for ID ${meta.id}: template=${templateCircle.base}, metadata=${metaBase}`);
                        return;
                    }
                    
                    // Verify coordinates (allow small floating point differences)
                    const xDiff = Math.abs(templateCircle.x - meta.x);
                    const yDiff = Math.abs(templateCircle.y - meta.y);
                    if (xDiff > 0.01 || yDiff > 0.01) {
                        warnings.push(`Coordinate mismatch for ID ${meta.id}: template=(${templateCircle.x},${templateCircle.y}), metadata=(${meta.x},${meta.y})`);
                    }
                    
                    verifiedCount++;
                });
                
                const result = {
                    valid: errors.length === 0,
                    verified: verifiedCount,
                    total: metadata.length,
                    errors: errors,
                    warnings: warnings
                };
                
                if (result.valid) {
                    console.log(`Metadata verification PASSED: ${verifiedCount}/${metadata.length} nodes verified`, warnings.length > 0 ? {warnings} : '');
                } else {
                    console.error('Metadata verification FAILED:', result);
                }
                
                return result;
            }
        };
        
        // Layer 2 & 3: METADATA-DRIVEN ARCHITECTURE with Verification
        // CRITICAL: Uses getElementById (O(1)) instead of querySelector (O(n)) for precise DOM access
        const RNAView = {
            metadata: RNA_METADATA || [],
            templateReference: null,  // Read-only SVG template reference
            registry: new Map(),  // id -> DOM Element direct reference (no selector ambiguity)
            buckets: { A: [], U: [], C: [], G: [] },  // Pre-indexed buckets: base -> [{dom, reactivity, id}]
            verificationResult: null,  // Store verification result
            
            init() {
                // Step 1: Extract and store SVG template reference (read-only)
                this.templateReference = SVGTemplateValidator.extractTemplateReference();
                
                // Step 2: Verify metadata against template
                if (this.metadata.length > 0 && this.templateReference) {
                    this.verificationResult = SVGTemplateValidator.verifyMetadata(this.metadata, this.templateReference);
                    
                    if (!this.verificationResult.valid) {
                        console.error('CRITICAL: Metadata verification failed. Visualization may be incorrect.');
                        // Continue anyway, but log errors
                    }
                }
                
                // Step 3: Clear previous state
                this.registry.clear();
                this.buckets = { A: [], U: [], C: [], G: [] };
                
                // Step 4: Build shadow map: physical ID -> DOM element reference
                // This happens ONCE at initialization, eliminating all future DOM queries
                this.metadata.forEach(meta => {
                    const normalizedBase = meta.base.toUpperCase();
                    // Normalize T to U (metadata should already have normalized bases, but ensure)
                    const base = normalizedBase === 'T' ? 'U' : normalizedBase;
                    
                    // ATOMIC DOM ACCESS: Use getElementById with physical ID (node_{id})
                    // This is O(1) and eliminates all selector ambiguity
                    const nodeId = `node_${meta.id}`;
                    const element = document.getElementById(nodeId);
                    
                    if (element) {
                        // Additional verification: cross-check with template reference
                        if (this.templateReference) {
                            const templateCircle = this.templateReference.circles.get(nodeId);
                            if (templateCircle) {
                                // Verify base matches template
                                if (templateCircle.base !== base) {
                                    console.error(`Base mismatch for ${nodeId}: metadata=${base}, template=${templateCircle.base}`);
                                }
                            }
                        }
                        
                        // Store direct DOM reference in registry (no future queries needed)
                        this.registry.set(meta.id, element);
                        
                        // Build bucket: pre-index by base type with direct DOM references
                        // This creates physical isolation - A bucket cannot contain G elements
                        if (this.buckets[base]) {
                            this.buckets[base].push({
                                dom: element,           // Direct DOM reference
                                reactivity: meta.reactivity,
                                id: meta.id,
                                base: base             // Store base for additional verification
                            });
                        }
                        
                        // Ensure data-base attribute is set correctly (for legacy compatibility)
                        element.setAttribute('data-base', base);
                    } else {
                        console.warn(`Element with ID ${nodeId} not found in DOM`);
                    }
                });
                
                // Step 5: Final bucket verification - ensure no cross-contamination
                const bucketVerification = this.verifyBucketIntegrity();
                if (!bucketVerification.valid) {
                    console.error('CRITICAL: Bucket integrity check failed:', bucketVerification.errors);
                }
                
                console.log('RNAView initialized with verification:', {
                    totalNodes: this.metadata.length,
                    registrySize: this.registry.size,
                    verificationPassed: this.verificationResult?.valid ?? 'N/A',
                    buckets: Object.keys(this.buckets).map(b => ({
                        base: b, 
                        count: this.buckets[b].length,
                        sampleIds: this.buckets[b].slice(0, 3).map(item => item.id)
                    }))
                });
            },
            
            // Verify bucket integrity - ensure no cross-base contamination
            verifyBucketIntegrity() {
                const errors = [];
                const baseCounts = { A: 0, U: 0, C: 0, G: 0 };
                
                Object.keys(this.buckets).forEach(base => {
                    this.buckets[base].forEach(item => {
                        // Count actual bases in each bucket
                        if (item.base === base) {
                            baseCounts[base]++;
                        } else {
                            errors.push(`Bucket ${base} contains element with base ${item.base} (ID: ${item.id})`);
                        }
                    });
                });
                
                return {
                    valid: errors.length === 0,
                    baseCounts: baseCounts,
                    errors: errors
                };
            },
            
            // Layer 4: ATOMIC UPDATE with Verification
            // Update colors for a specific base type using bucket-based approach
            // This function physically cannot affect other base types due to bucket isolation
            updateByBase(baseType) {
                // Input validation
                let normalizedBase = baseType.toUpperCase().trim();
                if (normalizedBase === 'T') {
                    normalizedBase = 'U';
                }
                
                if (!['A', 'U', 'C', 'G'].includes(normalizedBase)) {
                    console.warn('Invalid base type:', baseType);
                    return;
                }
                
                // Get bucket for this base type (physical isolation - no other bases can be in this array)
                const bucket = this.buckets[normalizedBase] || [];
                if (bucket.length === 0) {
                    console.warn(`No elements found for base type ${normalizedBase}`);
                    return; // No circles for this base type
                }
                
                // Verification: Ensure all items in bucket match the requested base type
                const mismatches = bucket.filter(item => item.base !== normalizedBase);
                if (mismatches.length > 0) {
                    console.error(`CRITICAL: Bucket ${normalizedBase} contains ${mismatches.length} mismatched elements:`, 
                        mismatches.map(m => ({id: m.id, base: m.base})));
                    // Filter out mismatches for safety
                    const validBucket = bucket.filter(item => item.base === normalizedBase);
                    if (validBucket.length === 0) {
                        console.error(`No valid elements in bucket ${normalizedBase} after filtering`);
                        return;
                    }
                    // Use filtered bucket
                    return this.updateBucket(validBucket, normalizedBase);
                }
                
                // All items verified, proceed with update
                this.updateBucket(bucket, normalizedBase);
            },
            
            // Internal method to update a verified bucket
            updateBucket(bucket, normalizedBase) {
            const baseEnabled = {
                'A': document.getElementById('baseEnableA')?.checked ?? true,
                'T': document.getElementById('baseEnableT')?.checked ?? true,
                    'U': document.getElementById('baseEnableT')?.checked ?? true,  // U uses T checkbox
                'C': document.getElementById('baseEnableC')?.checked ?? true,
                'G': document.getElementById('baseEnableG')?.checked ?? true
            };
            
            const strokeWidth = parseFloat(document.getElementById('strokeWidth')?.value || 2);
            
                // ATOMIC UPDATE: Only iterate over this base's bucket
                // Since buckets are physically separated and verified, this loop CANNOT access other base types
                let updatedCount = 0;
                bucket.forEach(item => {
                    // Final verification: ensure item base matches
                    if (item.base !== normalizedBase) {
                        console.error(`Skipping mismatched element: ID=${item.id}, expected=${normalizedBase}, actual=${item.base}`);
                        return;
                    }
                    
                    const element = item.dom;  // Direct DOM reference (no querySelector)
                    const reactivity = item.reactivity;
                    
                    // Additional verification: check element still exists and has correct ID
                    if (!element || !element.parentNode) {
                        console.warn(`Element for ID ${item.id} no longer exists in DOM`);
                        return;
                    }
                    
                    const elementId = element.getAttribute('id');
                    if (elementId !== `node_${item.id}`) {
                        console.error(`ID mismatch: expected node_${item.id}, found ${elementId}`);
                        return;
                    }
                    
                    // Update color based on reactivity and current thresholds
                    const enabled = baseEnabled[normalizedBase];
                    
                    if (!enabled) {
                        element.style.display = 'none';
                    } else {
                        element.style.display = '';
                        const color = determineColor(reactivity, normalizedBase);
                        const currentStyle = element.getAttribute('style') || '';
                        let newStyle = currentStyle.replace(/fill:[^;]+;?/g, '')
                                                   .replace(/stroke:[^;]+;?/g, '')
                                                   .replace(/stroke-width:[^;]+;?/g, '');
                        if (circleFilled) {
                            newStyle += `fill: ${color}; stroke: ${color}; stroke-width: ${strokeWidth};`;
                        } else {
                            newStyle += `fill: none; stroke: ${color}; stroke-width: ${strokeWidth};`;
                        }
                        element.setAttribute('style', newStyle.trim());
                        updatedCount++;
                    }
                });
                
                // Log update result for debugging
                if (updatedCount !== bucket.length) {
                    console.warn(`Updated ${updatedCount}/${bucket.length} elements for base ${normalizedBase}`);
                }
            }
        };
        
        // Legacy function for backward compatibility - now uses metadata-driven approach
        function updateColorsForBase(baseType) {
            RNAView.updateByBase(baseType);
        }
        
        // Helper function to update a single circle's color
        // baseType is the normalized base type (A, T, C, G) that this circle should have
        function updateCircleColor(circle, baseType, baseEnabled, strokeWidth) {
            // Get position and reactivity from circle attributes
                let position = null;
                let reactivity = null;
                
                if (circle.hasAttribute('data-position')) {
                    position = parseInt(circle.getAttribute('data-position'));
                    if (circle.hasAttribute('data-reactivity')) {
                        reactivity = parseFloat(circle.getAttribute('data-reactivity'));
                    }
                } else {
                    // Try to get from parent's title element
                    const parent = circle.parentElement;
                    if (parent) {
                        const title = parent.querySelector('title');
                        if (title) {
                            const match = title.textContent.match(/(\d+)/);
                            if (match) {
                                position = parseInt(match[1]);
                                circle.setAttribute('data-position', position);
                            }
                        }
                    }
                }
                
                // Get reactivity from data if available
                if (position && reactivityData.positions && reactivityData.positions[position]) {
                    const data = reactivityData.positions[position];
                    reactivity = data.reactivity;
                    // Set data attributes if not already set
                    if (!circle.hasAttribute('data-base')) {
                        circle.setAttribute('data-base', data.base);
                        circle.setAttribute('data-reactivity', data.reactivity);
                    }
                }
                
                if (!position || reactivity === null) {
                    // If no position or no reactivity data, skip (preserve original appearance)
                    return;
                }
                
            // Verify that the circle's base matches the expected baseType
            // This is a safety check since we already filtered by data-base attribute
            let circleBase = circle.getAttribute('data-base');
            if (circleBase) {
                circleBase = circleBase.toUpperCase().trim();
                if (circleBase === 'U') {
                    circleBase = 'T';
                }
                // If the circle's base doesn't match expected baseType, skip it
                // (This shouldn't happen if selector is correct, but safety check)
                if (circleBase !== baseType) {
                    return;
                }
            }
            
            // Use the baseType parameter (already normalized) for color determination
            const enabled = baseEnabled[baseType];
                
                if (!enabled) {
                    // Hide circles for disabled bases
                    circle.style.display = 'none';
                } else {
                    circle.style.display = '';
                // Update color based on reactivity using current color ranges for this specific base
                const color = determineColor(reactivity, baseType);
                    const currentStyle = circle.getAttribute('style') || '';
                    // Update style attribute, preserving other styles
                    let newStyle = currentStyle.replace(/fill:[^;]+;?/g, '').replace(/stroke:[^;]+;?/g, '').replace(/stroke-width:[^;]+;?/g, '');
                    if (circleFilled) {
                        newStyle += `fill: ${color}; stroke: ${color}; stroke-width: ${strokeWidth};`;
                    } else {
                        newStyle += `fill: none; stroke: ${color}; stroke-width: ${strokeWidth};`;
                    }
                    circle.setAttribute('style', newStyle.trim());
                }
        }
        
        // Apply colors and visibility based on current settings (for all bases)
        // METADATA-DRIVEN: Uses RNAView for precise, atomic updates
        function applyColorAndVisibility() {
            // Use metadata-driven approach if available
            if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                // Ensure registry is initialized
                if (RNAView.registry.size === 0) {
                    RNAView.init();
                }
                // Update each base type using metadata-driven approach
                ['A', 'U', 'C', 'G'].forEach(base => {
                    RNAView.updateByBase(base);
                });
            } else {
                // Fallback to legacy approach if metadata not available
                const baseEnabled = {
                    'A': document.getElementById('baseEnableA')?.checked ?? true,
                    'T': document.getElementById('baseEnableT')?.checked ?? true,
                    'C': document.getElementById('baseEnableC')?.checked ?? true,
                    'G': document.getElementById('baseEnableG')?.checked ?? true
                };
                
                const strokeWidth = parseFloat(document.getElementById('strokeWidth')?.value || 2);
                
                // Update each base type separately to ensure correct matching
                ['A', 'T', 'C', 'G'].forEach(base => {
                    updateColorsForBase(base);
                });
            }
        }
        
        // Update circle style (filled vs hollow)
        function updateCircleStyle() {
            const style = document.getElementById('circleStyle').value;
            circleFilled = (style === 'filled');
            
            // Apply style change by calling applyColorAndVisibility
            applyColorAndVisibility();
        }
        
        // Update font size and color - only for base letters (A, T, C, G, U)
        function updateFontStyle() {
            const fontSizeInput = document.getElementById('fontSize');
            const fontColorInput = document.getElementById('fontColor');
            
            if (!fontSizeInput || !fontColorInput) return;
            
            const fontSize = fontSizeInput.value + 'px';
            const fontColor = getColorValue('fontColor');
            
            // Update color preview
            const fontColorPreview = document.getElementById('fontColorPreview');
            if (fontColorPreview) {
                fontColorPreview.style.backgroundColor = fontColor;
            }
            
            // Convert hex to rgb/rgba if needed
            let finalColor = fontColor;
            if (fontColor.startsWith('#')) {
                const r = parseInt(fontColor.slice(1, 3), 16);
                const g = parseInt(fontColor.slice(3, 5), 16);
                const b = parseInt(fontColor.slice(5, 7), 16);
                finalColor = `rgb(${r}, ${g}, ${b})`;
            }
            
            // Only update text elements that are base letters (inside <g><title>position</title><text>base</text></g>)
            const allTexts = document.querySelectorAll('#svgContainer svg text');
            allTexts.forEach(text => {
                // Check if this text is a base letter (has parent <g> with <title> containing position number)
                const parent = text.parentElement;
                if (parent && parent.tagName === 'g') {
                    const title = parent.querySelector('title');
                    if (title) {
                        const titleText = title.textContent || '';
                        // Check if title contains position number and text is a single base letter
                        const textContent = text.textContent.trim();
                        if (/^\d+/.test(titleText) && /^[ATCGU]$/i.test(textContent)) {
                            // This is a base letter, update it
                            text.setAttribute('font-size', fontSize);
                            text.style.fontSize = fontSize;
                            text.setAttribute('fill', finalColor);
                            text.style.fill = finalColor;
                        }
                    }
                }
            });
        }
        
        // Update circle radius
        function updateCircleRadius() {
            const radiusInput = document.getElementById('circleRadius');
            if (!radiusInput) return;
            
            const radius = parseFloat(radiusInput.value);
            const radiusValueEl = document.getElementById('circleRadiusValue');
            if (radiusValueEl) {
                radiusValueEl.textContent = radius.toFixed(1);
            }
            
            // Update all circles
            const circles = document.querySelectorAll('#svgContainer svg circle');
            circles.forEach(circle => {
                // Only update circles that have position data (not legend circles)
                if (circle.hasAttribute('data-position') || 
                    (circle.parentElement && circle.parentElement.querySelector('title'))) {
                    circle.setAttribute('r', radius);
                }
            });
        }
        
        // Update legend position and direction - control SVG internal legend
        function updateLegend() {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            const position = document.getElementById('legendPosition')?.value || 'bottom';
            const direction = document.getElementById('legendDirection')?.value || 'horizontal';
            const fontSize = parseFloat(document.getElementById('legendFontSize')?.value || 12);
            const itemWidth = parseFloat(document.getElementById('legendItemWidth')?.value || 30);
            const itemHeight = parseFloat(document.getElementById('legendItemHeight')?.value || 15);
            const textAlign = document.getElementById('legendTextAlign')?.value || 'middle';
            const textOffsetX = parseFloat(document.getElementById('legendTextOffsetX')?.value || 0);
            const textOffsetY = parseFloat(document.getElementById('legendTextOffsetY')?.value || 0);
            const rowSpacing = parseFloat(document.getElementById('legendRowSpacing')?.value || 40);
            const colSpacing = parseFloat(document.getElementById('legendColSpacing')?.value || 50);
            
            // Find legend groups - use the new structured format
            // Look for groups with class "legend-group" or id starting with "legend-"
            let legends = Array.from(svg.querySelectorAll('g.legend-group, g[id^="legend-"], g[id^="color-bar-"]'));
            
            // Fallback: try to find by data-base attribute or other patterns
            if (legends.length === 0) {
                legends = Array.from(svg.querySelectorAll('g[data-base], g[class*="legend" i], g[id*="legend" i]'));
            }
            
            // If still no legend found, try to find by text content (old format)
            if (legends.length === 0) {
                const allGroups = svg.querySelectorAll('g');
                allGroups.forEach(g => {
                    const text = g.textContent || '';
                    if (text.toLowerCase().includes('legend') || 
                        text.toLowerCase().includes('color') ||
                        (g.querySelectorAll('rect, circle').length > 0 && g.querySelectorAll('text').length > 0)) {
                        legends.push(g);
                    }
                });
            }
            
            // If still no legend, look for groups with multiple colored rectangles/circles and text
            if (legends.length === 0) {
                const allGroups = svg.querySelectorAll('g');
                allGroups.forEach(g => {
                    const shapes = g.querySelectorAll('rect, circle');
                    const texts = g.querySelectorAll('text');
                    if (shapes.length >= 2 && texts.length >= 2) {
                        // Likely a legend
                        legends.push(g);
                    }
                });
            }
            
            // Group legends by base type if possible, or process individually
            legends.forEach((legend, legendIdx) => {
                try {
                    // Check if this is the new structured format
                    const itemsGroup = legend.querySelector('.legend-items');
                    const items = itemsGroup ? Array.from(itemsGroup.children) : Array.from(legend.children);
                    
                    const svgWidth = svg.viewBox.baseVal.width || svg.getBBox().width;
                    const svgHeight = svg.viewBox.baseVal.height || svg.getBBox().height;
                    
                    // Get legend bounding box
                    let legendBBox;
                    try {
                        legendBBox = legend.getBBox();
                    } catch (e) {
                        // Estimate from items
                        let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
                        items.forEach(item => {
                            try {
                                const bbox = item.getBBox();
                                minX = Math.min(minX, bbox.x);
                                minY = Math.min(minY, bbox.y);
                                maxX = Math.max(maxX, bbox.x + bbox.width);
                                maxY = Math.max(maxY, bbox.y + bbox.height);
                            } catch (e) {}
                        });
                        legendBBox = {x: minX, y: minY, width: maxX - minX, height: maxY - minY};
                    }
                    
                    // Calculate layout for multiple legends (one per base)
                    // Horizontal direction: each legend (base) on its own row (stacked vertically)
                    // Vertical direction: each legend (base) in its own column (side by side) or stacked in one column with spacing
                    
                    // First rearrange items within legend, then calculate dimensions
                    // We'll calculate actual dimensions after items are rearranged
                    
                    // Use spacing from controls (already set above)
                    
                    // We'll calculate actual legend dimensions after items rearrangement
                    // Use configured item dimensions
                    const estimatedItemWidth = itemWidth; // Width of each color segment (horizontal)
                    const estimatedItemHeight = itemHeight; // Height of each color segment (vertical)
                    const itemsCount = items.length;
                    
                    let legendWidth, legendHeight;
                    if (direction === 'horizontal') {
                        // Horizontal: legend width is sum of all item widths, height is item height
                        legendWidth = itemsCount * estimatedItemWidth;
                        legendHeight = estimatedItemHeight;
                    } else {
                        // Vertical: legend width is max item width, height is sum of all item heights
                        legendWidth = estimatedItemWidth;
                        legendHeight = itemsCount * estimatedItemHeight;
                    }
                    
                    // Calculate total space needed for all legends
                    let totalLegendsWidth, totalLegendsHeight;
                    if (direction === 'horizontal') {
                        // Horizontal: legends stacked vertically (each base on its own row)
                        totalLegendsWidth = legendWidth; // All legends same width
                        totalLegendsHeight = legends.length * legendHeight + (legends.length - 1) * rowSpacing;
                    } else {
                        // Vertical: legends side by side (each base in its own column)
                        totalLegendsWidth = legends.length * legendWidth + (legends.length - 1) * colSpacing;
                        totalLegendsHeight = legendHeight; // All legends same height
                    }
                    
                    let baseX = 0, baseY = 0;
                    
                    // Apply position - calculate base position for each legend group
                    if (position === 'top') {
                        if (direction === 'horizontal') {
                            // Horizontal: each legend on its own row, centered horizontally
                            baseX = svgWidth / 2 - legendWidth / 2;
                            baseY = 20 + legendIdx * (legendHeight + rowSpacing);
                        } else {
                            // Vertical: each legend in its own column, centered vertically
                            baseX = svgWidth / 2 - totalLegendsWidth / 2 + legendIdx * (legendWidth + colSpacing);
                            baseY = 20;
                        }
                    } else if (position === 'bottom') {
                        if (direction === 'horizontal') {
                            // Horizontal: each legend on its own row, centered horizontally
                            baseX = svgWidth / 2 - legendWidth / 2;
                            baseY = svgHeight - totalLegendsHeight - 20 + legendIdx * (legendHeight + rowSpacing);
                        } else {
                            // Vertical: each legend in its own column, centered vertically
                            baseX = svgWidth / 2 - totalLegendsWidth / 2 + legendIdx * (legendWidth + colSpacing);
                            baseY = svgHeight - legendHeight - 20;
                        }
                    } else if (position === 'left') {
                        if (direction === 'horizontal') {
                            // Horizontal: each legend on its own row, aligned to left
                            baseX = 20;
                            baseY = svgHeight / 2 - totalLegendsHeight / 2 + legendIdx * (legendHeight + rowSpacing);
                        } else {
                            // Vertical: each legend in its own column, aligned to left, centered vertically
                            baseX = 20 + legendIdx * (legendWidth + colSpacing);
                            baseY = svgHeight / 2 - legendHeight / 2;
                        }
                    } else if (position === 'right') {
                        if (direction === 'horizontal') {
                            // Horizontal: each legend on its own row, aligned to right
                            baseX = svgWidth - legendWidth - 20;
                            baseY = svgHeight / 2 - totalLegendsHeight / 2 + legendIdx * (legendHeight + rowSpacing);
                        } else {
                            // Vertical: each legend in its own column, aligned to right, centered vertically
                            baseX = svgWidth - totalLegendsWidth - 20 + legendIdx * (legendWidth + colSpacing);
                            baseY = svgHeight / 2 - legendHeight / 2;
                        }
                    }
                    
                    // Apply font size to legend text elements
                    const legendTexts = legend.querySelectorAll('text');
                    legendTexts.forEach(text => {
                        text.setAttribute('font-size', fontSize);
                        text.style.fontSize = fontSize + 'px';
                    });
                    
                    // First, save original item positions and dimensions BEFORE any transforms
                    // We need to get the original positions by temporarily removing transforms
                    const originalTransforms = new Map();
                    const originalPositions = new Map();
                    const originalDimensions = new Map();
                    
                    // Store current transforms and clear them temporarily to get original positions
                    items.forEach((item, idx) => {
                        // Store current transform
                        const currentTransform = item.getAttribute('transform') || '';
                        originalTransforms.set(item, currentTransform);
                        
                        // Temporarily remove transform to get original bbox
                        item.removeAttribute('transform');
                        
                        try {
                            const bbox = item.getBBox();
                            originalPositions.set(item, { x: bbox.x, y: bbox.y });
                            originalDimensions.set(item, { width: bbox.width, height: bbox.height });
                        } catch (e) {
                            // If getBBox fails, estimate from structure
                            const rect = item.querySelector('rect');
                            const circle = item.querySelector('circle');
                            const text = item.querySelector('text');
                            let width = 20, height = 20;
                            if (rect) {
                                width = parseFloat(rect.getAttribute('width') || 20);
                                height = parseFloat(rect.getAttribute('height') || 20);
                            } else if (circle) {
                                const r = parseFloat(circle.getAttribute('r') || 10);
                                width = height = r * 2;
                            }
                            originalPositions.set(item, { x: idx * 30, y: 0 });
                            originalDimensions.set(item, { width, height });
                        }
                    });
                    
                    // Rearrange elements within each item based on direction
                    // New structured format: each item is a <g class="legend-item"> containing rect and text
                    items.forEach((item) => {
                        // Skip if not a group element
                        if (item.tagName !== 'g') return;
                        
                        const rect = item.querySelector('rect.legend-rect, rect');
                        const text = item.querySelector('text.legend-text, text');
                        
                        if (!rect || !text) return;
                        
                        // Use configured dimensions
                        const finalRectWidth = itemWidth;
                        const finalRectHeight = itemHeight;
                        const textSize = parseFloat(text.getAttribute('font-size') || fontSize);
                        
                        // Update rect dimensions to match configured values
                        rect.setAttribute('width', finalRectWidth.toString());
                        rect.setAttribute('height', finalRectHeight.toString());
                        
                        if (direction === 'horizontal') {
                            // Horizontal layout: rect should be wide (horizontal bar), text above
                            // Position rect: x relative to item, y below text
                            rect.setAttribute('x', '0');
                            rect.setAttribute('y', (textSize + 2).toString()); // Below text with 2px spacing
                            rect.removeAttribute('transform');
                            
                            // Position text above rect based on alignment, with offset
                            let textX = finalRectWidth / 2; // Default to middle
                            if (textAlign === 'start') {
                                textX = 0;
                            } else if (textAlign === 'end') {
                                textX = finalRectWidth;
                            }
                            textX += textOffsetX; // Apply X offset
                            const textY = textSize + textOffsetY; // Apply Y offset
                            text.setAttribute('x', textX.toString());
                            text.setAttribute('y', textY.toString());
                            text.setAttribute('text-anchor', textAlign);
                            text.setAttribute('font-size', fontSize.toString());
                            text.removeAttribute('transform');
                        } else {
                            // Vertical layout: rect should be tall (vertical bar), text below
                            // Position rect at top
                            rect.setAttribute('x', '0');
                            rect.setAttribute('y', '0');
                            rect.removeAttribute('transform');
                            
                            // Position text below rect based on alignment, with offset
                            let textX = finalRectWidth / 2; // Default to middle
                            if (textAlign === 'start') {
                                textX = 0;
                            } else if (textAlign === 'end') {
                                textX = finalRectWidth;
                            }
                            textX += textOffsetX; // Apply X offset
                            const textY = finalRectHeight + textSize + 2 + textOffsetY; // Below rect with spacing and Y offset
                            text.setAttribute('x', textX.toString());
                            text.setAttribute('y', textY.toString());
                            text.setAttribute('text-anchor', textAlign);
                            text.setAttribute('font-size', fontSize.toString());
                            text.removeAttribute('transform');
                        }
                    });
                    
                    // Now rearrange items based on direction
                    // For horizontal: align all items to the same Y coordinate (baseline)
                    // For vertical: align all items to the same X coordinate (baseline)
                    
                    // Recalculate item dimensions after internal rearrangement
                    // IMPORTANT: Read rect dimensions AFTER they've been adjusted by direction
                    const itemDimensions = new Map();
                    items.forEach((item) => {
                        const rect = item.querySelector('rect.legend-rect, rect');
                        const text = item.querySelector('text.legend-text, text');
                        
                        if (rect) {
                            // Read rect dimensions AFTER adjustment (they should already be correct orientation)
                            const rectWidth = parseFloat(rect.getAttribute('width') || itemWidth);
                            const rectHeight = parseFloat(rect.getAttribute('height') || itemHeight);
                            
                            if (direction === 'horizontal') {
                                // Horizontal: width is rect width (for continuous bar)
                                // Height includes text above + rect
                                const textSize = text ? parseFloat(text.getAttribute('font-size') || fontSize) : fontSize;
                                itemDimensions.set(item, { 
                                    width: rectWidth,  // Only rect width for continuous horizontal bar
                                    height: textSize + 2 + rectHeight  // text + spacing + rect
                                });
                            } else {
                                // Vertical: height is rect height (for continuous bar)
                                // Width includes text and rect
                                const textBBox = text ? text.getBBox() : { width: 0 };
                                const textWidth = textBBox.width || rectWidth;
                                itemDimensions.set(item, { 
                                    width: Math.max(rectWidth, textWidth),  // Max of rect and text width
                                    height: rectHeight  // Only rect height for continuous vertical bar
                                });
                            }
                        } else {
                            // Fallback: use bbox
                            try {
                                const itemTransform = item.getAttribute('transform') || '';
                                item.removeAttribute('transform');
                                const bbox = item.getBBox();
                                itemDimensions.set(item, { width: bbox.width, height: bbox.height });
                                item.setAttribute('transform', itemTransform);
                            } catch (e) {
                                const dim = originalDimensions.get(item);
                                if (dim) {
                                    itemDimensions.set(item, dim);
                                }
                            }
                        }
                    });
                    
                    if (direction === 'horizontal') {
                        // Arrange items horizontally in a row, all aligned to Y=0
                        // Items form a continuous horizontal bar (rects connect directly)
                        let currentX = 0;
                        items.forEach((item) => {
                            const dim = itemDimensions.get(item);
                            if (dim) {
                                // All items start at Y=0, creating a continuous horizontal bar
                                item.setAttribute('transform', `translate(${currentX}, 0)`);
                                // Move to next position - use rect width only for continuous bar
                                currentX += dim.width; // No spacing - rects connect directly
                            }
                        });
                    } else {
                        // Arrange items vertically in a column, all aligned to X=0
                        // Items form a continuous vertical bar (rects connect directly)
                        let currentY = 0;
                        items.forEach((item) => {
                            const dim = itemDimensions.get(item);
                            if (dim) {
                                // All items start at X=0, creating a continuous vertical bar
                                item.setAttribute('transform', `translate(0, ${currentY})`);
                                // Move to next position - use rect height only for continuous bar
                                currentY += dim.height; // No spacing - rects connect directly
                            }
                        });
                    }
                    
                    // Position base label relative to items based on direction
                    const baseLabel = legend.querySelector('.legend-base-label, text[class*="base-label" i]');
                    if (baseLabel && items.length > 0) {
                        // Get items group bounding box to position base label relative to it
                        const itemsGroup = legend.querySelector('.legend-items');
                        let itemsBBox;
                        try {
                            if (itemsGroup) {
                                itemsBBox = itemsGroup.getBBox();
                            } else {
                                // Fallback: calculate from items
                                let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
                                items.forEach(item => {
                                    try {
                                        const bbox = item.getBBox();
                                        minX = Math.min(minX, bbox.x);
                                        minY = Math.min(minY, bbox.y);
                                        maxX = Math.max(maxX, bbox.x + bbox.width);
                                        maxY = Math.max(maxY, bbox.y + bbox.height);
                                    } catch (e) {}
                                });
                                itemsBBox = {x: minX, y: minY, width: maxX - minX, height: maxY - minY};
                            }
                        } catch (e) {
                            // Estimate from item dimensions
                            const firstItemDim = itemDimensions.get(items[0]);
                            if (firstItemDim) {
                                if (direction === 'horizontal') {
                                    itemsBBox = {x: 0, y: 0, width: items.length * itemWidth, height: itemHeight + fontSize + 2};
                                } else {
                                    itemsBBox = {x: 0, y: 0, width: itemWidth, height: items.length * itemHeight + fontSize + 2};
                                }
                            } else {
                                itemsBBox = {x: 0, y: 0, width: 100, height: 20};
                            }
                        }
                        
                        // Position base label based on direction
                        if (direction === 'horizontal') {
                            // Horizontal: base label on the left, vertically centered with items
                            const labelY = itemsBBox.y + itemsBBox.height / 2;
                            baseLabel.setAttribute('x', (itemsBBox.x - 20).toString()); // 20px spacing to the left
                            baseLabel.setAttribute('y', labelY.toString());
                            baseLabel.setAttribute('text-anchor', 'end'); // Right-align text (since it's on the left)
                            baseLabel.setAttribute('dominant-baseline', 'middle'); // Vertically center
                        } else {
                            // Vertical: base label on top, horizontally centered with items
                            const labelX = itemsBBox.x + itemsBBox.width / 2;
                            baseLabel.setAttribute('x', labelX.toString());
                            baseLabel.setAttribute('y', (itemsBBox.y - 5).toString()); // 5px spacing above
                            baseLabel.setAttribute('text-anchor', 'middle'); // Center horizontally
                            baseLabel.setAttribute('dominant-baseline', 'baseline'); // Align to baseline
                        }
                        // Apply font size to base label
                        baseLabel.setAttribute('font-size', fontSize.toString());
                    }
                    
                    // Now get the updated bounding box after rearranging items
                    let updatedLegendBBox;
                    try {
                        // Temporarily remove group transform to get accurate bbox
                        const groupTransform = legend.getAttribute('transform') || '';
                        legend.removeAttribute('transform');
                        updatedLegendBBox = legend.getBBox();
                        legend.setAttribute('transform', groupTransform);
                    } catch (e) {
                        // Fallback: calculate from items and base label
                        let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
                        items.forEach(item => {
                            try {
                                const bbox = item.getBBox();
                                minX = Math.min(minX, bbox.x);
                                minY = Math.min(minY, bbox.y);
                                maxX = Math.max(maxX, bbox.x + bbox.width);
                                maxY = Math.max(maxY, bbox.y + bbox.height);
                            } catch (e) {}
                        });
                        // Include base label in bbox calculation
                        if (baseLabel) {
                            try {
                                const labelBBox = baseLabel.getBBox();
                                minX = Math.min(minX, labelBBox.x);
                                minY = Math.min(minY, labelBBox.y);
                                maxX = Math.max(maxX, labelBBox.x + labelBBox.width);
                                maxY = Math.max(maxY, labelBBox.y + labelBBox.height);
                            } catch (e) {}
                        }
                        updatedLegendBBox = {x: minX, y: minY, width: maxX - minX, height: maxY - minY};
                    }
                    
                    // Finally, apply the group transform to position the legend
                    legend.setAttribute('transform', `translate(${baseX - updatedLegendBBox.x}, ${baseY - updatedLegendBBox.y})`);
                } catch (e) {
                    console.warn('Error updating legend:', e);
                }
            });
        }
        
        // Initialize legend dragging functionality - drag all legends together as a group
        function initLegendDragging() {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            let isDraggingLegend = false;
            let dragStart = { x: 0, y: 0 };
            let legendStartTransforms = new Map(); // Store initial transforms for all legends
            
            // Find all legend groups
            const legends = Array.from(svg.querySelectorAll('g.legend-group, g[id^="legend-"]'));
            
            if (legends.length === 0) return;
            
            // Create a wrapper group or use the first legend as reference
            // Make all legends draggable by adding cursor style and event listeners
            legends.forEach(legend => {
                // Make legend draggable by adding cursor style
                legend.style.cursor = 'move';
                
                // Mouse down - start dragging all legends together
                legend.addEventListener('mousedown', (e) => {
                    // Don't drag if clicking on text or other interactive elements
                    if (e.target.tagName === 'text' || e.target.closest('text')) {
                        return;
                    }
                    
                    // Don't drag if clicking on control panel or other UI elements
                    if (e.target.closest('.control-panel')) {
                        return;
                    }
                    
                    isDraggingLegend = true;
                    
                    // Store initial transforms for all legends
                    legendStartTransforms.clear();
                    legends.forEach(l => {
                        const currentTransform = l.getAttribute('transform') || '';
                        const match = currentTransform.match(/translate\(([^,]+),\s*([^)]+)\)/);
                        if (match) {
                            legendStartTransforms.set(l, {
                                x: parseFloat(match[1]) || 0,
                                y: parseFloat(match[2]) || 0
                            });
                        } else {
                            legendStartTransforms.set(l, { x: 0, y: 0 });
                        }
                    });
                    
                    // Get mouse position in SVG coordinates
                    const svgPoint = svg.createSVGPoint();
                    svgPoint.x = e.clientX;
                    svgPoint.y = e.clientY;
                    const ctm = svg.getScreenCTM();
                    if (ctm) {
                        const svgCTM = ctm.inverse();
                        const point = svgPoint.matrixTransform(svgCTM);
                        dragStart.x = point.x;
                        dragStart.y = point.y;
                    } else {
                        dragStart.x = e.clientX;
                        dragStart.y = e.clientY;
                    }
                    
                    e.preventDefault();
                    e.stopPropagation();
                });
            });
            
            // Mouse move - update all legend positions together
            svg.addEventListener('mousemove', (e) => {
                if (!isDraggingLegend || legendStartTransforms.size === 0) return;
                
                // Get mouse position in SVG coordinates
                const svgPoint = svg.createSVGPoint();
                svgPoint.x = e.clientX;
                svgPoint.y = e.clientY;
                const ctm = svg.getScreenCTM();
                if (ctm) {
                    const svgCTM = ctm.inverse();
                    const point = svgPoint.matrixTransform(svgCTM);
                    
                    // Calculate delta
                    const deltaX = point.x - dragStart.x;
                    const deltaY = point.y - dragStart.y;
                    
                    // Update transform for all legends
                    legendStartTransforms.forEach((startTransform, legend) => {
                        const newX = startTransform.x + deltaX;
                        const newY = startTransform.y + deltaY;
                        legend.setAttribute('transform', `translate(${newX}, ${newY})`);
                    });
                }
                
                e.preventDefault();
            });
            
            // Mouse up - stop dragging
            document.addEventListener('mouseup', () => {
                if (isDraggingLegend) {
                    isDraggingLegend = false;
                    legendStartTransforms.clear();
                }
            });
            
            // Also handle mouse leave to stop dragging
            svg.addEventListener('mouseleave', () => {
                if (isDraggingLegend) {
                    isDraggingLegend = false;
                    legendStartTransforms.clear();
                }
            });
        }
        
        // Minimap auto-hide timer
        let minimapHideTimer = null;
        let minimapActivityTimer = null;
        const MINIMAP_HIDE_DELAY = 3000; // 3 seconds
        
        // Reset minimap auto-hide timer (called when user activity stops)
        function resetMinimapHideTimer() {
            // Clear existing timer
            if (minimapHideTimer) {
                clearTimeout(minimapHideTimer);
                minimapHideTimer = null;
            }
            
            // Clear activity debounce timer
            if (minimapActivityTimer) {
                clearTimeout(minimapActivityTimer);
            }
            
            const minimap = document.getElementById('minimap');
            if (!minimap) return;
            
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            const viewBox = svg.viewBox.baseVal;
            
            // Only set timer if minimap should be visible (zoomed in)
            if (viewBox.width < initialViewBox.width * 0.9) {
                // Use debounce: wait a short time after last activity before starting hide timer
                minimapActivityTimer = setTimeout(() => {
                    minimap.style.display = 'block';
                    // Start hide timer only after activity has stopped
                    minimapHideTimer = setTimeout(() => {
                        minimap.style.display = 'none';
                        minimapHideTimer = null;
                    }, MINIMAP_HIDE_DELAY);
                }, 100); // 100ms debounce to avoid flickering
            }
        }
        
        // Update minimap (only updates display, doesn't reset timer)
        function updateMinimap() {
            const svg = document.querySelector('#svgContainer svg');
            const minimap = document.getElementById('minimap');
            const minimapSvg = document.getElementById('minimapSvg');
            const viewport = document.getElementById('minimapViewport');
            
            if (!svg || !minimap || !minimapSvg || !viewport) return;
            
            const viewBox = svg.viewBox.baseVal;
            const svgRect = svg.getBBox();
            
            // Show minimap only when zoomed
            if (viewBox.width < initialViewBox.width * 0.9) {
                minimap.style.display = 'block';
                
                // Calculate viewport position in minimap
                const scaleX = 100 / initialViewBox.width;
                const scaleY = 100 / initialViewBox.height;
                
                const x = ((viewBox.x - initialViewBox.x) * scaleX);
                const y = ((viewBox.y - initialViewBox.y) * scaleY);
                const width = (viewBox.width * scaleX);
                const height = (viewBox.height * scaleY);
                
                viewport.setAttribute('x', Math.max(0, x));
                viewport.setAttribute('y', Math.max(0, y));
                viewport.setAttribute('width', Math.min(100, width));
                viewport.setAttribute('height', Math.min(100, height));
                
                // Don't reset timer here - only update display
                // Timer will be reset by user activity handlers
            } else {
                minimap.style.display = 'none';
                if (minimapHideTimer) {
                    clearTimeout(minimapHideTimer);
                    minimapHideTimer = null;
                }
                if (minimapActivityTimer) {
                    clearTimeout(minimapActivityTimer);
                    minimapActivityTimer = null;
                }
            }
        }
        
        // Minimap click and drag handler
        function setupMinimap() {
            const minimapSvg = document.getElementById('minimapSvg');
            const viewport = document.getElementById('minimapViewport');
            const svg = document.querySelector('#svgContainer svg');
            
            if (!minimapSvg || !svg || !viewport) return;
            
            let isDraggingViewport = false;
            let dragStart = { x: 0, y: 0 };
            let viewportStart = { x: 0, y: 0 };
            
            // Make viewport draggable - use capture phase to ensure event is handled
            viewport.addEventListener('mousedown', (e) => {
                e.preventDefault();
                e.stopPropagation();
                isDraggingViewport = true;
                
                const rect = minimapSvg.getBoundingClientRect();
                dragStart.x = e.clientX;
                dragStart.y = e.clientY;
                
                // Get current viewport position in minimap coordinates (0-100)
                viewportStart.x = parseFloat(viewport.getAttribute('x')) || 0;
                viewportStart.y = parseFloat(viewport.getAttribute('y')) || 0;
                
                // Show minimap and cancel any hide timer
                const minimap = document.getElementById('minimap');
                if (minimap) {
                    minimap.style.display = 'block';
                }
                if (minimapHideTimer) {
                    clearTimeout(minimapHideTimer);
                    minimapHideTimer = null;
                }
                if (minimapActivityTimer) {
                    clearTimeout(minimapActivityTimer);
                    minimapActivityTimer = null;
                }
            }, true); // Use capture phase
            
            // Handle viewport dragging - use capture phase
            const handleViewportDrag = (e) => {
                if (!isDraggingViewport) return;
                
                e.preventDefault();
                
                const rect = minimapSvg.getBoundingClientRect();
                const minimapWidth = rect.width;
                const minimapHeight = rect.height;
                
                // Calculate mouse movement in minimap coordinates (0-100)
                const deltaX = ((e.clientX - dragStart.x) / minimapWidth) * 100;
                const deltaY = ((e.clientY - dragStart.y) / minimapHeight) * 100;
                
                // Calculate new viewport position
                let newX = viewportStart.x + deltaX;
                let newY = viewportStart.y + deltaY;
                
                // Get viewport dimensions
                const viewportWidth = parseFloat(viewport.getAttribute('width')) || 100;
                const viewportHeight = parseFloat(viewport.getAttribute('height')) || 100;
                
                // Constrain to minimap bounds (0-100)
                newX = Math.max(0, Math.min(100 - viewportWidth, newX));
                newY = Math.max(0, Math.min(100 - viewportHeight, newY));
                
                // Update viewport position
                viewport.setAttribute('x', newX);
                viewport.setAttribute('y', newY);
                
                // Update main SVG viewBox based on viewport position
                const scaleX = initialViewBox.width / 100;
                const scaleY = initialViewBox.height / 100;
                
                const viewBox = svg.viewBox.baseVal;
                viewBox.x = initialViewBox.x + (newX * scaleX);
                viewBox.y = initialViewBox.y + (newY * scaleY);
                
                // Constrain viewBox to bounds
                viewBox.x = Math.max(initialViewBox.x, Math.min(initialViewBox.x + initialViewBox.width - viewBox.width, viewBox.x));
                viewBox.y = Math.max(initialViewBox.y, Math.min(initialViewBox.y + initialViewBox.height - viewBox.height, viewBox.y));
                
                updateZoomDisplayFromViewBox();
                
                // Cancel any hide timer while dragging
                if (minimapHideTimer) {
                    clearTimeout(minimapHideTimer);
                    minimapHideTimer = null;
                }
                if (minimapActivityTimer) {
                    clearTimeout(minimapActivityTimer);
                    minimapActivityTimer = null;
                }
            };
            
            document.addEventListener('mousemove', handleViewportDrag, true);
            
            // Stop dragging on mouse up
            const handleViewportDragEnd = () => {
                if (isDraggingViewport) {
                    isDraggingViewport = false;
                    // Start hide timer after dragging ends (no activity for 3 seconds)
                    resetMinimapHideTimer();
                }
            };
            
            document.addEventListener('mouseup', handleViewportDragEnd, true);
            document.addEventListener('mouseleave', handleViewportDragEnd, true);
            
            // Click on minimap background to jump to that position
            minimapSvg.addEventListener('click', (e) => {
                // Don't handle click if clicking on viewport (handled by drag)
                if (e.target === viewport) return;
                
                const rect = minimapSvg.getBoundingClientRect();
                const x = ((e.clientX - rect.left) / rect.width) * 100;
                const y = ((e.clientY - rect.top) / rect.height) * 100;
                
                const viewBox = svg.viewBox.baseVal;
                const newX = (x / 100) * initialViewBox.width + initialViewBox.x - viewBox.width / 2;
                const newY = (y / 100) * initialViewBox.height + initialViewBox.y - viewBox.height / 2;
                
                viewBox.x = Math.max(initialViewBox.x, Math.min(initialViewBox.x + initialViewBox.width - viewBox.width, newX));
                viewBox.y = Math.max(initialViewBox.y, Math.min(initialViewBox.y + initialViewBox.height - viewBox.height, newY));
                
                updateMinimap();
                updateZoomDisplayFromViewBox();
                
                // Start hide timer after click (no activity for 3 seconds)
                resetMinimapHideTimer();
            });
        }
        
        
        // Reset view function
        function resetView() {
            const svg = document.querySelector('#svgContainer svg');
            if (svg && initialViewBox) {
                const viewBox = svg.viewBox.baseVal;
                viewBox.x = initialViewBox.x;
                viewBox.y = initialViewBox.y;
                viewBox.width = initialViewBox.width;
                viewBox.height = initialViewBox.height;
            }
            document.getElementById('circleStyle').value = circleFilled ? 'filled' : 'hollow';
            
            // Reset color ranges to initial values
            document.getElementById('rangeA1').value = 0.3;
            document.getElementById('rangeA2').value = 0.65;
            document.getElementById('rangeT1').value = 0.2;
            document.getElementById('rangeT2').value = 0.55;
            document.getElementById('rangeC1').value = 0.3;
            document.getElementById('rangeC2').value = 0.6;
            document.getElementById('rangeG1').value = 0.2;
            document.getElementById('rangeG2').value = 0.5;
            
            // Reset base enables
            document.getElementById('baseEnableA').checked = true;
            document.getElementById('baseEnableT').checked = true;
            document.getElementById('baseEnableC').checked = true;
            document.getElementById('baseEnableG').checked = true;
            
            updateColorRanges();
            updateRangeDisplay();
            updateMinimap();
            updateZoomDisplayFromViewBox();
        }
        
        // Update range display values
        function updateRangeDisplay() {
            ['A', 'T', 'C', 'G'].forEach(base => {
                const range1El = document.getElementById(`range${base}1`);
                const range2El = document.getElementById(`range${base}2`);
                const range1ValueEl = document.getElementById(`range${base}1Value`);
                const range2ValueEl = document.getElementById(`range${base}2Value`);
                
                if (range1El && range1ValueEl) {
                    const value = parseFloat(range1El.value);
                    range1ValueEl.textContent = value.toFixed(2);
                    console.log(`Updated range${base}1 display:`, value);
                } else {
                    console.warn(`Missing elements for range${base}1`);
                }
                if (range2El && range2ValueEl) {
                    const value = parseFloat(range2El.value);
                    range2ValueEl.textContent = value.toFixed(2);
                    console.log(`Updated range${base}2 display:`, value);
                } else {
                    console.warn(`Missing elements for range${base}2`);
                }
            });
        }
        
        // Export function with format support
        function exportImage() {
            const format = document.getElementById('exportFormat').value;
            const fileName = document.getElementById('exportFileName').value || 'rna_structure_interactive';
            const svg = document.querySelector('#svgContainer svg');
            
            if (!svg) return;
            
            if (format === 'svg') {
                const serializer = new XMLSerializer();
                const svgString = serializer.serializeToString(svg);
                const blob = new Blob([svgString], { type: 'image/svg+xml' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = fileName + '.svg';
                a.click();
                URL.revokeObjectURL(url);
            } else if (format === 'png') {
                const resolution = parseFloat(document.getElementById('pngResolution').value);
                const dpiEl = document.getElementById('pngDPI');
                const dpi = dpiEl ? parseFloat(dpiEl.value) : 72;
                const svgData = new XMLSerializer().serializeToString(svg);
                const canvas = document.createElement('canvas');
                const ctx = canvas.getContext('2d');
                const img = new Image();
                
                const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
                const url = URL.createObjectURL(svgBlob);
                
                img.onload = function() {
                    // Calculate dimensions based on DPI
                    const scale = dpi / 72; // 72 is standard screen DPI
                    canvas.width = img.width * resolution * scale;
                    canvas.height = img.height * resolution * scale;
                    ctx.scale(resolution * scale, resolution * scale);
                    ctx.drawImage(img, 0, 0);
                    canvas.toBlob(function(blob) {
                        const a = document.createElement('a');
                        a.href = URL.createObjectURL(blob);
                        a.download = fileName + '.png';
                        a.click();
                        URL.revokeObjectURL(a.href);
                    }, 'image/png');
                    URL.revokeObjectURL(url);
                };
                img.src = url;
            } else if (format === 'pdf') {
                // PDF export requires jsPDF library - simplified version
                alert('PDF export requires additional library. Please use SVG or PNG format.');
            }
        }
        
        // Legacy exportSVG function for backward compatibility
        function exportSVG() {
            document.getElementById('exportFormat').value = 'svg';
            exportImage();
        }
        
        // Undo/Redo functionality
        let undoStack = [];
        let redoStack = [];
        let maxUndoSteps = 50;
        
        function saveStateForUndo() {
            const state = {
                circles: Array.from(document.querySelectorAll('#svgContainer svg circle')).map(c => ({
                    position: c.getAttribute('data-position'),
                    r: c.getAttribute('r'),
                    style: c.getAttribute('style')
                })),
                texts: Array.from(document.querySelectorAll('#svgContainer svg text')).map(t => ({
                    text: t.textContent,
                    fontSize: t.getAttribute('font-size'),
                    fill: t.getAttribute('fill')
                }))
            };
            undoStack.push(state);
            if (undoStack.length > maxUndoSteps) {
                undoStack.shift();
            }
            redoStack = []; // Clear redo stack when new action is performed
            updateUndoRedoButtons();
        }
        
        function undoAction() {
            if (undoStack.length === 0) return;
            const currentState = {
                circles: Array.from(document.querySelectorAll('#svgContainer svg circle')).map(c => ({
                    position: c.getAttribute('data-position'),
                    r: c.getAttribute('r'),
                    style: c.getAttribute('style')
                })),
                texts: Array.from(document.querySelectorAll('#svgContainer svg text')).map(t => ({
                    text: t.textContent,
                    fontSize: t.getAttribute('font-size'),
                    fill: t.getAttribute('fill')
                }))
            };
            redoStack.push(currentState);
            const previousState = undoStack.pop();
            restoreState(previousState);
            updateUndoRedoButtons();
        }
        
        function redoAction() {
            if (redoStack.length === 0) return;
            const currentState = {
                circles: Array.from(document.querySelectorAll('#svgContainer svg circle')).map(c => ({
                    position: c.getAttribute('data-position'),
                    r: c.getAttribute('r'),
                    style: c.getAttribute('style')
                })),
                texts: Array.from(document.querySelectorAll('#svgContainer svg text')).map(t => ({
                    text: t.textContent,
                    fontSize: t.getAttribute('font-size'),
                    fill: t.getAttribute('fill')
                }))
            };
            undoStack.push(currentState);
            const nextState = redoStack.pop();
            restoreState(nextState);
            updateUndoRedoButtons();
        }
        
        function restoreState(state) {
            if (!state) return;
            
            // Restore circles
            if (state.circles) {
                state.circles.forEach(circleData => {
                    if (circleData.position) {
                        const circle = document.querySelector(`circle[data-position="${circleData.position}"]`);
                        if (circle) {
                            if (circleData.r !== undefined && circleData.r !== null) {
                                circle.setAttribute('r', circleData.r);
                            }
                            if (circleData.style) {
                                circle.setAttribute('style', circleData.style);
                            }
                        }
                    }
                });
            }
            
            // Restore texts
            if (state.texts) {
                state.texts.forEach(textData => {
                    const texts = Array.from(document.querySelectorAll('#svgContainer svg text'));
                    const matchingText = texts.find(t => {
                        const textContent = t.textContent.trim();
                        return textContent === textData.text || 
                               (textData.text && textContent === textData.text.trim());
                    });
                    if (matchingText) {
                        if (textData.fontSize) {
                            matchingText.setAttribute('font-size', textData.fontSize);
                            matchingText.style.fontSize = textData.fontSize;
                        }
                        if (textData.fill) {
                            matchingText.setAttribute('fill', textData.fill);
                            matchingText.style.fill = textData.fill;
                        }
                    }
                });
                
                // Update font size input to match restored state
                const firstText = state.texts.find(t => t.fontSize);
                if (firstText && firstText.fontSize) {
                    const fontSizeValue = parseFloat(firstText.fontSize.replace('px', ''));
                    if (!isNaN(fontSizeValue)) {
                        const fontSizeInput = document.getElementById('fontSize');
                        if (fontSizeInput) {
                            fontSizeInput.value = fontSizeValue;
                            const fontSizeValueEl = document.getElementById('fontSizeValue');
                            if (fontSizeValueEl) {
                                fontSizeValueEl.textContent = fontSizeValue.toFixed(1);
                            }
                        }
                    }
                }
                
                // Update font color input if available
                const firstTextWithColor = state.texts.find(t => t.fill);
                if (firstTextWithColor && firstTextWithColor.fill) {
                    const fontColorInput = document.getElementById('fontColor');
                    if (fontColorInput) {
                        // Convert rgb/rgba to hex if needed
                        let colorValue = firstTextWithColor.fill;
                        if (colorValue.startsWith('rgb')) {
                            const rgb = colorValue.match(/\d+/g);
                            if (rgb && rgb.length >= 3) {
                                const r = parseInt(rgb[0]).toString(16).padStart(2, '0');
                                const g = parseInt(rgb[1]).toString(16).padStart(2, '0');
                                const b = parseInt(rgb[2]).toString(16).padStart(2, '0');
                                colorValue = '#' + r + g + b;
                            }
                        }
                        fontColorInput.value = colorValue;
                        updateColorPreview('fontColor', 'fontColorPreview');
                        if (typeof $ !== 'undefined' && $(fontColorInput).spectrum) {
                            $(fontColorInput).spectrum('set', colorValue);
                        }
                    }
                }
            }
        }
        
        function updateUndoRedoButtons() {
            const undoBtn = document.getElementById('undoBtn');
            const redoBtn = document.getElementById('redoBtn');
            if (undoBtn) undoBtn.disabled = undoStack.length === 0;
            if (redoBtn) redoBtn.disabled = redoStack.length === 0;
        }
        
        // Update zoom display
        function updateZoomDisplay() {
            updateZoomDisplayFromViewBox();
        }
        
        // Event listeners
        document.getElementById('circleStyle').addEventListener('change', updateCircleStyle);
        document.getElementById('strokeWidth').addEventListener('input', () => {
            document.getElementById('strokeWidthValue').textContent = 
                parseFloat(document.getElementById('strokeWidth').value).toFixed(1);
            applyColorAndVisibility();
        });
        // Track if we're currently adjusting to avoid saving state on every input event
        let isAdjustingFontSize = false;
        let isAdjustingCircleRadius = false;
        
        const fontSizeEl = document.getElementById('fontSize');
        if (fontSizeEl) {
            fontSizeEl.addEventListener('input', () => {
                const fontSizeValueEl = document.getElementById('fontSizeValue');
                if (fontSizeValueEl) {
                    fontSizeValueEl.textContent = fontSizeEl.value;
                }
                updateFontStyle();
                isAdjustingFontSize = true;
            });
            fontSizeEl.addEventListener('mousedown', () => {
                // Save state when starting to adjust
                saveStateForUndo();
            });
            fontSizeEl.addEventListener('mouseup', () => {
                // Save state when finished adjusting
                if (isAdjustingFontSize) {
                    saveStateForUndo();
                    isAdjustingFontSize = false;
                }
            });
            fontSizeEl.addEventListener('change', () => {
                // Also save on change event (for keyboard input)
                if (isAdjustingFontSize) {
                    saveStateForUndo();
                    isAdjustingFontSize = false;
                }
            });
        }
        
        const circleRadiusEl = document.getElementById('circleRadius');
        if (circleRadiusEl) {
            circleRadiusEl.addEventListener('input', () => {
                updateCircleRadius();
                isAdjustingCircleRadius = true;
            });
            circleRadiusEl.addEventListener('mousedown', () => {
                // Save state when starting to adjust
                saveStateForUndo();
            });
            circleRadiusEl.addEventListener('mouseup', () => {
                // Save state when finished adjusting
                if (isAdjustingCircleRadius) {
                    saveStateForUndo();
                    isAdjustingCircleRadius = false;
                }
            });
            circleRadiusEl.addEventListener('change', () => {
                // Also save on change event (for keyboard input)
                if (isAdjustingCircleRadius) {
                    saveStateForUndo();
                    isAdjustingCircleRadius = false;
                }
            });
        }
        
        const fontColorEl = document.getElementById('fontColor');
        if (fontColorEl) {
            // Handle direct text input changes
            fontColorEl.addEventListener('input', (e) => {
                const colorValue = e.target.value;
                // Validate hex color format
                if (/^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/.test(colorValue)) {
                    updateFontStyle();
                    updateColorPreview('fontColor', 'fontColorPreview');
                    // Update Spectrum if initialized
                    if (typeof $ !== 'undefined' && $(fontColorEl).spectrum) {
                        try {
                            $(fontColorEl).spectrum('set', colorValue);
                        } catch (err) {
                            console.warn('Failed to update Spectrum:', err);
                        }
                    }
                    saveState();
                }
            });
            fontColorEl.addEventListener('blur', () => {
                // On blur, validate and correct the color value
                const colorValue = fontColorEl.value;
                if (!/^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/.test(colorValue)) {
                    // Invalid format, restore previous valid value
                    const prevValue = getColorValue('fontColor') || '#000000';
                    fontColorEl.value = prevValue;
                    updateFontStyle();
                    updateColorPreview('fontColor', 'fontColorPreview');
                }
            });
            // Listen for Spectrum color picker changes
            if (typeof $ !== 'undefined') {
                // Wait for Spectrum to be initialized, then bind the change event
                setTimeout(() => {
                    $(fontColorEl).off('change.fontColor'); // Remove any existing handlers
                    $(fontColorEl).on('change.fontColor', () => {
                        updateFontStyle();
                        updateColorPreview('fontColor', 'fontColorPreview');
                        saveState();
                    });
                }, 200);
            }
            
            // Make color preview clickable to open color picker
            const fontColorPreview = document.getElementById('fontColorPreview');
            if (fontColorPreview) {
                fontColorPreview.addEventListener('click', (e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    if (typeof $ !== 'undefined' && $(fontColorEl).spectrum) {
                        $(fontColorEl).spectrum('show');
                    } else if (fontColorEl) {
                        fontColorEl.focus();
                        fontColorEl.select();
                    }
                });
            }
        }
        
        // Background color handling
        const backgroundColorEl = document.getElementById('backgroundColor');
        if (backgroundColorEl) {
            // Handle direct text input changes
            backgroundColorEl.addEventListener('input', (e) => {
                const colorValue = e.target.value;
                // Validate hex color format
                if (/^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/.test(colorValue)) {
                    updateBackgroundColor();
                    updateColorPreview('backgroundColor', 'backgroundColorPreview');
                    // Update Spectrum if initialized
                    if (typeof $ !== 'undefined' && $(backgroundColorEl).spectrum) {
                        try {
                            $(backgroundColorEl).spectrum('set', colorValue);
                        } catch (err) {
                            console.warn('Failed to update Spectrum:', err);
                        }
                    }
                    saveState();
                }
            });
            backgroundColorEl.addEventListener('blur', () => {
                // On blur, validate and correct the color value
                const colorValue = backgroundColorEl.value;
                if (!/^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/.test(colorValue)) {
                    // Invalid format, restore previous valid value
                    const prevValue = getColorValue('backgroundColor') || '#ffffff';
                    backgroundColorEl.value = prevValue;
                    updateBackgroundColor();
                    updateColorPreview('backgroundColor', 'backgroundColorPreview');
                }
            });
            // Listen for Spectrum color picker changes (use 'change' event, not 'change.spectrum')
            if (typeof $ !== 'undefined') {
                // Wait for Spectrum to be initialized, then bind the change event
                setTimeout(() => {
                    $(backgroundColorEl).off('change.backgroundColor'); // Remove any existing handlers
                    $(backgroundColorEl).on('change.backgroundColor', () => {
                        updateBackgroundColor();
                        updateColorPreview('backgroundColor', 'backgroundColorPreview');
                        saveState();
                    });
                }, 200);
            }
            
            // Make color preview clickable to open color picker
            const backgroundColorPreview = document.getElementById('backgroundColorPreview');
            if (backgroundColorPreview) {
                backgroundColorPreview.addEventListener('click', () => {
                    if (typeof $ !== 'undefined' && $(backgroundColorEl).spectrum) {
                        $(backgroundColorEl).spectrum('show');
                    } else {
                        backgroundColorEl.click();
                    }
                });
            }
        }
        
        // Update color previews for all color pickers
        function updateColorPreview(inputId, previewId) {
            const input = document.getElementById(inputId);
            const preview = document.getElementById(previewId);
            if (input && preview) {
                try {
                    // First try to use input.value directly (most reliable, especially during initialization)
                    if (input.value) {
                        // Check if it's a valid hex color
                        if (/^#[0-9A-Fa-f]{6}$/i.test(input.value)) {
                            preview.style.backgroundColor = input.value;
                            return;
                        }
                        // If it's a valid color format (even if not hex), try to use it
                        if (/^#[0-9A-Fa-f]{3,8}$/i.test(input.value) || /^rgb|rgba|hsl|hsla/i.test(input.value)) {
                            preview.style.backgroundColor = input.value;
                            return;
                        }
                    }
                    
                    // Only try Spectrum if jQuery is available and input has a value
                    if (typeof $ !== 'undefined' && input.value) {
                        try {
                            // Check if Spectrum is initialized
                            const spectrumInstance = $(input).data('spectrum');
                            if (spectrumInstance) {
                                const color = getColorValue(inputId);
                                if (color && color !== '#000000') {
                                    preview.style.backgroundColor = color;
                                    return;
                                }
                            }
                        } catch (e) {
                            // Ignore Spectrum errors, fall back to input value
                        }
                    }
                    
                    // Fallback to input value or default
                    if (input.value) {
                        preview.style.backgroundColor = input.value;
                    } else {
                        preview.style.backgroundColor = '#000000';
                    }
                } catch (e) {
                    console.warn('Error updating color preview:', e, 'inputId:', inputId);
                    // Fallback to input value or default
                    if (input && input.value) {
                        preview.style.backgroundColor = input.value;
                    } else {
                        preview.style.backgroundColor = '#000000';
                    }
                }
            }
        }
        
        // Initialize color previews
        function initColorPreviews() {
            ['A', 'T', 'C', 'G'].forEach(base => {
                for (let i = 1; i <= 3; i++) {
                    const inputId = `color${base}${i}`;
                    const previewId = `${inputId}Preview`;
                    const input = document.getElementById(inputId);
                    const preview = document.getElementById(previewId);
                    if (input && preview) {
                        preview.style.backgroundColor = input.value;
                        input.addEventListener('input', () => updateColorPreview(inputId, previewId));
                        if (typeof $ !== 'undefined') {
                            $(input).on('change', () => updateColorPreview(inputId, previewId));
                        }
                    }
                }
            });
            // Unified colors
            ['unifiedColor1', 'unifiedColor2', 'unifiedColor3'].forEach(id => {
                const previewId = `${id}Preview`;
                const input = document.getElementById(id);
                const preview = document.getElementById(previewId);
                if (input && preview) {
                    preview.style.backgroundColor = input.value;
                    input.addEventListener('input', () => updateColorPreview(id, previewId));
                    if (typeof $ !== 'undefined') {
                        $(input).on('change', () => updateColorPreview(id, previewId));
                    }
                }
            });
            // Update font color preview - use direct value during initialization to avoid Spectrum issues
            const fontColorInput = document.getElementById('fontColor');
            const fontColorPreviewEl = document.getElementById('fontColorPreview');
            if (fontColorInput && fontColorPreviewEl) {
                fontColorPreviewEl.style.backgroundColor = fontColorInput.value || '#000000';
            }
            // Update background color preview
            const backgroundColorInput = document.getElementById('backgroundColor');
            const backgroundColorPreviewEl = document.getElementById('backgroundColorPreview');
            if (backgroundColorInput && backgroundColorPreviewEl) {
                backgroundColorPreviewEl.style.backgroundColor = backgroundColorInput.value || '#ffffff';
            }
        }
        
        // Update background color
        function updateBackgroundColor() {
            const backgroundColorInput = document.getElementById('backgroundColor');
            const svgContainer = document.getElementById('svgContainer');
            if (!backgroundColorInput || !svgContainer) return;
            
            const color = getColorValue('backgroundColor') || '#ffffff';
            svgContainer.style.backgroundColor = color;
            
            // Update preview
            const backgroundColorPreview = document.getElementById('backgroundColorPreview');
            if (backgroundColorPreview) {
                backgroundColorPreview.style.backgroundColor = color;
            }
        }
        
        document.getElementById('legendPosition').addEventListener('change', updateLegend);
        document.getElementById('legendDirection').addEventListener('change', updateLegend);
        document.getElementById('legendTextAlign').addEventListener('change', updateLegend);
        
        // Text offset controls
        const legendTextOffsetXEl = document.getElementById('legendTextOffsetX');
        if (legendTextOffsetXEl) {
            legendTextOffsetXEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendTextOffsetXValue');
                if (valueEl) {
                    valueEl.textContent = legendTextOffsetXEl.value;
                }
                updateLegend();
            });
        }
        
        const legendTextOffsetYEl = document.getElementById('legendTextOffsetY');
        if (legendTextOffsetYEl) {
            legendTextOffsetYEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendTextOffsetYValue');
                if (valueEl) {
                    valueEl.textContent = legendTextOffsetYEl.value;
                }
                updateLegend();
            });
        }
        
        // Spacing controls
        const legendRowSpacingEl = document.getElementById('legendRowSpacing');
        if (legendRowSpacingEl) {
            legendRowSpacingEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendRowSpacingValue');
                if (valueEl) {
                    valueEl.textContent = legendRowSpacingEl.value;
                }
                updateLegend();
            });
        }
        
        const legendColSpacingEl = document.getElementById('legendColSpacing');
        if (legendColSpacingEl) {
            legendColSpacingEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendColSpacingValue');
                if (valueEl) {
                    valueEl.textContent = legendColSpacingEl.value;
                }
                updateLegend();
            });
        }
        
        const legendFontSizeEl = document.getElementById('legendFontSize');
        if (legendFontSizeEl) {
            legendFontSizeEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendFontSizeValue');
                if (valueEl) {
                    valueEl.textContent = legendFontSizeEl.value;
                }
                updateLegend();
            });
            legendFontSizeEl.addEventListener('change', updateLegend);
        }
        const legendItemWidthEl = document.getElementById('legendItemWidth');
        if (legendItemWidthEl) {
            legendItemWidthEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendItemWidthValue');
                if (valueEl) {
                    valueEl.textContent = legendItemWidthEl.value;
                }
                updateLegend();
            });
            legendItemWidthEl.addEventListener('change', updateLegend);
        }
        const legendItemHeightEl = document.getElementById('legendItemHeight');
        if (legendItemHeightEl) {
            legendItemHeightEl.addEventListener('input', () => {
                const valueEl = document.getElementById('legendItemHeightValue');
                if (valueEl) {
                    valueEl.textContent = legendItemHeightEl.value;
                }
                updateLegend();
            });
            legendItemHeightEl.addEventListener('change', updateLegend);
        }
        
        // Export format change handler
        const exportFormatEl = document.getElementById('exportFormat');
        if (exportFormatEl) {
            exportFormatEl.addEventListener('change', () => {
                const format = exportFormatEl.value;
                const pngResolutionRow = document.getElementById('pngResolutionRow');
                const pngDPIRow = document.getElementById('pngDPIRow');
                if (format === 'png') {
                    if (pngResolutionRow) pngResolutionRow.style.display = '';
                    if (pngDPIRow) pngDPIRow.style.display = '';
                } else {
                    if (pngResolutionRow) pngResolutionRow.style.display = 'none';
                    if (pngDPIRow) pngDPIRow.style.display = 'none';
                }
            });
        }
        
        // Color range controls
        ['A', 'T', 'C', 'G'].forEach(base => {
            const range1El = document.getElementById(`range${base}1`);
            const range2El = document.getElementById(`range${base}2`);
            const baseEnableEl = document.getElementById(`baseEnable${base}`);
            
            if (range1El) {
                range1El.addEventListener('input', () => {
                    updateRangeDisplay();
                    // CRITICAL FIX: Only update color range for this specific base (not all bases)
                    updateColorRangeForBase(base);
                    updateLegendValues(base); // CRITICAL FIX: Only update legend for this specific base
                    // CRITICAL FIX: Only update colors for this specific base
                    updateColorsForBase(base);
                    // OPTIMIZATION: Only update distribution preview for this specific base
                    updateDistributionPreviews(base);
                });
            }
            if (range2El) {
                range2El.addEventListener('input', () => {
                    updateRangeDisplay();
                    // CRITICAL FIX: Only update color range for this specific base (not all bases)
                    updateColorRangeForBase(base);
                    updateLegendValues(base); // CRITICAL FIX: Only update legend for this specific base
                    // CRITICAL FIX: Only update colors for this specific base
                    updateColorsForBase(base);
                    // OPTIMIZATION: Only update distribution preview for this specific base
                    updateDistributionPreviews(base);
                });
            }
            if (baseEnableEl) {
                baseEnableEl.addEventListener('change', applyColorAndVisibility);
            }
            
            // Update distribution preview when color changes
            for (let i = 1; i <= 3; i++) {
                const colorEl = document.getElementById(`color${base}${i}`);
                if (colorEl) {
                    if (typeof $ !== 'undefined') {
                        $(colorEl).on('change', () => {
                            // CRITICAL FIX: Only update color range for this specific base (not all bases)
                            updateColorRangeForBase(base);
                            updateLegendValues(base); // CRITICAL FIX: Only update legend for this specific base
                            // CRITICAL FIX: Only update colors for this specific base
                            updateColorsForBase(base);
                            // OPTIMIZATION: Only update distribution preview for this specific base
                            updateDistributionPreviews(base);
                        });
                    } else {
                        colorEl.addEventListener('change', () => {
                            // CRITICAL FIX: Only update color range for this specific base (not all bases)
                            updateColorRangeForBase(base);
                            updateLegendValues(base); // CRITICAL FIX: Only update legend for this specific base
                            // CRITICAL FIX: Only update colors for this specific base
                            updateColorsForBase(base);
                            // OPTIMIZATION: Only update distribution preview for this specific base
                            updateDistributionPreviews(base);
                        });
                    }
                }
            }
        });
        
        // Unified controls event listeners - with verification
        const unifiedRange1El = document.getElementById('unifiedRange1');
        const unifiedRange2El = document.getElementById('unifiedRange2');
        
        if (unifiedRange1El) {
            unifiedRange1El.addEventListener('input', () => {
                const valueEl = document.getElementById('unifiedRange1Value');
                if (valueEl) {
                    valueEl.textContent = parseFloat(unifiedRange1El.value).toFixed(2);
                }
                updateColorRanges();
                // CRITICAL: Update all bases using verified metadata-driven approach
                if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                    ['A', 'U', 'C', 'G'].forEach(base => {
                        RNAView.updateByBase(base);
                    });
                } else {
                    applyColorAndVisibility();
                }
            });
        }
        if (unifiedRange2El) {
            unifiedRange2El.addEventListener('input', () => {
                const valueEl = document.getElementById('unifiedRange2Value');
                if (valueEl) {
                    valueEl.textContent = parseFloat(unifiedRange2El.value).toFixed(2);
                }
                updateColorRanges();
                // CRITICAL: Update all bases using verified metadata-driven approach
                if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                    ['A', 'U', 'C', 'G'].forEach(base => {
                        RNAView.updateByBase(base);
                    });
                } else {
                    applyColorAndVisibility();
                }
            });
        }
        ['unifiedColor1', 'unifiedColor2', 'unifiedColor3'].forEach(id => {
            const el = document.getElementById(id);
            if (el) {
                if (typeof $ !== 'undefined') {
                    $(el).on('change', () => {
                        updateColorRanges();
                        updateDistributionPreviews(); // Update distribution preview when color changes
                        // CRITICAL: Update all bases using verified metadata-driven approach
                        if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                            ['A', 'U', 'C', 'G'].forEach(base => {
                                RNAView.updateByBase(base);
                            });
                        } else {
                            applyColorAndVisibility();
                        }
                    });
                } else {
                    el.addEventListener('change', () => {
                        updateColorRanges();
                        updateDistributionPreviews(); // Update distribution preview when color changes
                        // CRITICAL: Update all bases using verified metadata-driven approach
                        if (RNAView && RNAView.metadata && RNAView.metadata.length > 0) {
                            ['A', 'U', 'C', 'G'].forEach(base => {
                                RNAView.updateByBase(base);
                    });
                        } else {
                            applyColorAndVisibility();
                        }
                    });
                }
            }
        });
        
        const unifiedColorRangesEl = document.getElementById('unifiedColorRanges');
        if (unifiedColorRangesEl) {
            unifiedColorRangesEl.addEventListener('change', () => {
                // toggleUnifiedMode() already handles updateColorRanges(), applyColorAndVisibility(), 
                // updateDistributionPreviews(), and updateLegendValues()
                toggleUnifiedMode();
            });
        }
        
        // Initialize distribution previews after page load
        setTimeout(() => {
            updateDistributionPreviews();
        }, 500);
        
        // Update distribution previews on window resize
        let resizeTimeout;
        window.addEventListener('resize', () => {
            clearTimeout(resizeTimeout);
            resizeTimeout = setTimeout(() => {
                updateDistributionPreviews();
            }, 200);
        });
        
        // Keyboard shortcuts
        document.addEventListener('keydown', (e) => {
            // Don't trigger shortcuts when typing in inputs
            if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') {
                if (e.key === 'Escape') {
                    document.getElementById('infoPanel').classList.remove('visible');
                    document.getElementById('searchPanel').classList.remove('visible');
                }
                return;
            }
            
            switch(e.key) {
                case '+':
                case '=':
                    e.preventDefault();
                    zoomIn();
                    break;
                case '-':
                case '_':
                    e.preventDefault();
                    zoomOut();
                    break;
                case '0':
                    e.preventDefault();
                    resetZoom();
                    break;
                case 'f':
                case 'F':
                    e.preventDefault();
                    resetZoom();
                    break;
                case 'r':
                case 'R':
                    e.preventDefault();
                    resetView();
                    break;
                case 'e':
                case 'E':
                    e.preventDefault();
                    exportImage();
                    break;
                case 'z':
                case 'Z':
                    if (e.ctrlKey || e.metaKey) {
                        e.preventDefault();
                        if (e.shiftKey) {
                            redoAction();
                        } else {
                            undoAction();
                        }
                    }
                    break;
                case 'y':
                case 'Y':
                    if (e.ctrlKey || e.metaKey) {
                        e.preventDefault();
                        redoAction();
                    }
                    break;
                case 'i':
                case 'I':
                    e.preventDefault();
                    showInfoPanel();
                    break;
                case 's':
                case 'S':
                    e.preventDefault();
                    showSearchPanel();
                    break;
                case 'h':
                case 'H':
                    e.preventDefault();
                    showKeyboardHelp();
                    break;
                case 'Escape':
                    document.getElementById('infoPanel').classList.remove('visible');
                    document.getElementById('searchPanel').classList.remove('visible');
                    break;
            }
        });
        
        // Save state on changes
        function setupStateSaving() {
            const saveableInputs = [
                'zoomLevel', 'circleStyle', 'strokeWidth', 'fontSize', 'fontColor', 'backgroundColor', 'unifiedColorRanges'
            ];
            
            saveableInputs.forEach(id => {
                const el = document.getElementById(id);
                if (el) {
                    // Background color and font color are handled by Spectrum change callback in initColorPickers
                    // and by separate event listeners in initializeVisualizationFeatures
                    // So we only need to handle saveState here
                    if (id === 'backgroundColor' || id === 'fontColor') {
                        // These will be handled by the event listeners set up in initializeVisualizationFeatures
                        // Just ensure state saving is triggered
                    } else {
                        el.addEventListener('change', saveState);
                        el.addEventListener('input', () => {
                            clearTimeout(window.saveTimeout);
                            window.saveTimeout = setTimeout(saveState, 500);
                        });
                    }
                }
            });
            
            ['A', 'T', 'C', 'G'].forEach(base => {
                const enableEl = document.getElementById(`baseEnable${base}`);
                if (enableEl) enableEl.addEventListener('change', saveState);
                for (let i = 1; i <= 2; i++) {
                    const rangeEl = document.getElementById(`range${base}${i}`);
                    if (rangeEl) {
                        rangeEl.addEventListener('input', () => {
                            clearTimeout(window.saveTimeout);
                            window.saveTimeout = setTimeout(saveState, 500);
                        });
                    }
                }
                for (let i = 1; i <= 3; i++) {
                    const input = document.getElementById(`color${base}${i}`);
                    if (input) {
                        if (typeof $ !== 'undefined' && $(input).spectrum) {
                            $(input).on('change', saveState);
                        } else {
                            input.addEventListener('change', saveState);
                        }
                    }
                }
            });
        }
        
        // Initialize Spectrum color pickers after DOM is ready
        function initColorPickersDelayed() {
            setTimeout(() => {
                initColorPickers();
            }, 100);
        }
        
        // Detect original font size from SVG and set range
        function detectOriginalFontSize() {
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) return;
            
            const allTexts = Array.from(svg.querySelectorAll('text'));
            let minFontSize = Infinity;
            let maxFontSize = -Infinity;
            let baseLetterFontSize = null;
            
            allTexts.forEach(text => {
                const parent = text.parentElement;
                if (parent && parent.tagName === 'g') {
                    const title = parent.querySelector('title');
                    if (title) {
                        const titleText = title.textContent || '';
                        const textContent = text.textContent.trim();
                        if (/^\d+/.test(titleText) && /^[ATCGU]$/i.test(textContent)) {
                            // This is a base letter
                            // Get font size from attribute first, then style, then computed style
                            let fontSize = parseFloat(text.getAttribute('font-size'));
                            if (isNaN(fontSize)) {
                                fontSize = parseFloat(text.style.fontSize);
                            }
                            if (isNaN(fontSize)) {
                                const computed = window.getComputedStyle(text);
                                fontSize = parseFloat(computed.fontSize);
                            }
                            if (!isNaN(fontSize) && fontSize > 0) {
                                if (baseLetterFontSize === null) baseLetterFontSize = fontSize;
                                minFontSize = Math.min(minFontSize, fontSize);
                                maxFontSize = Math.max(maxFontSize, fontSize);
                            }
                        }
                    }
                }
            });
            
            // If we found base letter font sizes, update the range
            if (baseLetterFontSize !== null) {
                const fontSizeInput = document.getElementById('fontSize');
                if (fontSizeInput) {
                    // Set range to 0.1x to 5x of original size
                    // Allow smaller than original size (down to 0.1x)
                    const minSize = Math.max(0.1, baseLetterFontSize * 0.1);
                    const maxSize = baseLetterFontSize * 5;
                    fontSizeInput.min = minSize;
                    fontSizeInput.max = maxSize;
                    fontSizeInput.step = 0.1;
                    // Always set initial value to template font size (don't check for existing value)
                    fontSizeInput.value = baseLetterFontSize;
                    const fontSizeValueEl = document.getElementById('fontSizeValue');
                    if (fontSizeValueEl) {
                        fontSizeValueEl.textContent = baseLetterFontSize.toFixed(1);
                    }
                }
            }
        }
        
        // Wrap initialization in a function that can be called after SVG load
        function initializeVisualizationFeatures() {
            // Initialize - wait for SVG to be ready
            setTimeout(() => {
            // Check if SVG exists before initializing
            const svg = document.querySelector('#svgContainer svg');
            if (!svg) {
                console.warn('SVG not found, skipping initialization');
                return;
            }
            
            // METADATA-DRIVEN ARCHITECTURE: Initialize RNAView registry
            if (typeof RNAView !== 'undefined' && RNAView.metadata && RNAView.metadata.length > 0) {
                RNAView.init();
                console.log('RNAView initialized with', RNAView.metadata.length, 'base nodes');
            } else {
                console.warn('RNA_METADATA not available, using legacy approach');
            }
            
            // First detect original font size BEFORE loading saved state
            // This ensures we preserve the original SVG font size
            detectOriginalFontSize();
            
            // Load saved state (if exists)
            const hasSavedState = localStorage.getItem('modDetectorState');
            if (hasSavedState) {
                loadSavedState();
            }
            
            initColorPickersDelayed();
            // Initialize color previews - use direct input.value to avoid Spectrum initialization issues
            initColorPreviews(); // Initialize color previews
            setupStateSaving();
            
            // Update range display values first (before any other updates)
            updateRangeDisplay();
            
            // Update unified range display values
            const unifiedRange1 = document.getElementById('unifiedRange1');
            const unifiedRange2 = document.getElementById('unifiedRange2');
            if (unifiedRange1 && document.getElementById('unifiedRange1Value')) {
                document.getElementById('unifiedRange1Value').textContent = 
                    parseFloat(unifiedRange1.value).toFixed(2);
            }
            if (unifiedRange2 && document.getElementById('unifiedRange2Value')) {
                document.getElementById('unifiedRange2Value').textContent = 
                    parseFloat(unifiedRange2.value).toFixed(2);
            }
            
            updateZoomDisplayFromViewBox(); // Use viewBox-based display
            // Don't call updateColorRanges() on init - preserve original SVG colors
            // Only update when user changes settings
            
            // Always update background color and font style (these should always apply)
            // Use setTimeout to ensure color pickers are initialized first
            setTimeout(() => {
                updateBackgroundColor();
                updateFontStyle();
            }, 150);
            
            // Only update styles if user has saved state
            // This preserves original SVG appearance on first load or after replacing SVG
            if (hasSavedState) {
                // User has saved state, apply it
                setTimeout(() => {
                    updateCircleStyle();
                    updateFontStyle();
                    updateCircleRadius();
                    updateLegend();
                }, 150);
            }
            
            updateInfoPanel();
            toggleUnifiedMode(); // Initialize unified mode display
            
            // Initialize legend dragging
            initLegendDragging();
            
            // Initialize minimap dragging
            setupMinimap();
            
            // Initialize distribution previews after a short delay to ensure canvas is ready
            setTimeout(() => {
                updateRangeDisplay(); // Ensure display values are updated
                updateDistributionPreviews();
            }, 200);
            
            // Initialize tool mode
            setToolMode('select');
            
            // Initialize export format display
            const exportFormatEl = document.getElementById('exportFormat');
            if (exportFormatEl) {
                exportFormatEl.dispatchEvent(new Event('change'));
            }
            
            // Save initial state for undo
            saveStateForUndo();
            updateUndoRedoButtons(); // Initialize undo/redo buttons
            
            // Update info panel periodically (clear old interval if exists)
            if (window.infoPanelInterval) {
                clearInterval(window.infoPanelInterval);
            }
            window.infoPanelInterval = setInterval(updateInfoPanel, 2000);
            }, 100);
        }
        
        // Initialize if SVG is already present
        if (document.querySelector('#svgContainer svg') && document.querySelector('#svgContainer svg').innerHTML.trim().length > 0) {
            initializeVisualizationFeatures();
        }
        
        // Listen for SVG load event
        document.addEventListener('svgLoaded', initializeVisualizationFeatures);
    </script>
</body>
</html>"###);

    // Write HTML to file
    let mut output = File::create(output_file)?;
    output.write_all(html.as_bytes())?;
    output.flush()?;
    
    Ok(())
}

/// Generate interactive HTML visualization from SVG template and reactivity data
/// This function should use the aligned SVG (with circles already drawn) to extract reactivity data
pub fn generate_interactive_visualization(
    reactivity_file: &str,
    svg_file: &str, // This should be the generated SVG with circles, not the template
    output_file: &str,
    bases_filter: &str,
    signal_type: &str,
    strand_filter: &str,
    color_ranges: Option<BaseColorRanges>,
    initial_circle_filled: bool,
) -> Result<(), Box<dyn Error>> {
    // Read the generated SVG (with circles and data attributes)
    let svg_content = std::fs::read_to_string(svg_file)?;
    
    // Extract reactivity data from the SVG's circles (which already have aligned positions)
    let reactivity_data = extract_reactivity_from_svg(&svg_content)?;
    
    // Extract complete base metadata for metadata-driven architecture
    let base_metadata = extract_base_metadata_from_svg(&svg_content)?;
    
    // Use provided color ranges or default
    let color_ranges = color_ranges.unwrap_or_default();
    
    // Generate interactive HTML with metadata
    generate_interactive_html(&svg_content, &reactivity_data, output_file, signal_type, &color_ranges, initial_circle_filled, Some(&base_metadata))?;
    
    Ok(())
}

/// Extract reactivity data from SVG circles that have data-position and data-reactivity attributes
/// CRITICAL FIX: Parse each <circle> tag individually to avoid context window ambiguity
/// This prevents matching attributes from adjacent circles
fn extract_reactivity_from_svg(svg_content: &str) -> Result<Vec<SvgReactivityData>, Box<dyn Error>> {
    let mut data = Vec::new();
    
    // Match complete <circle> tags (including self-closing and with closing tag)
    // Pattern: <circle ... /> or <circle ...></circle>
    let circle_regex = Regex::new(r#"<circle[^>]*?/>|<circle([^>]*?)>(?:[^<]*</circle>)?"#)?;
    let position_regex = Regex::new(r#"data-position="(\d+)""#)?;
    let base_regex = Regex::new(r#"data-base="([ATCGU])""#)?;
    let reactivity_regex = Regex::new(r#"data-reactivity="([^"]+)""#)?;
    
    // Find all circle tags and parse each one individually
    for circle_cap in circle_regex.captures_iter(svg_content) {
        // Get the circle tag content (either self-closing or with attributes)
        let circle_tag = if let Some(attrs) = circle_cap.get(1) {
            attrs.as_str()
        } else {
            // Self-closing tag: extract attributes from the full match
            let full_match = circle_cap.get(0).unwrap().as_str();
            // Remove <circle and /> to get just attributes
            full_match.trim_start_matches("<circle").trim_end_matches("/>").trim()
        };
        
        // Extract attributes from this specific circle tag only
        if let Some(pos_cap) = position_regex.captures(circle_tag) {
            if let Ok(position) = pos_cap[1].parse::<u32>() {
                let base = base_regex.captures(circle_tag)
                    .and_then(|c| c[1].chars().next())
                    .unwrap_or('N');
                let reactivity = reactivity_regex.captures(circle_tag)
                    .and_then(|c| c[1].parse::<f64>().ok())
                    .unwrap_or(0.0);
                
                data.push(SvgReactivityData {
                    position,
                    base: normalize_base(base),
                    reactivity,
                });
            }
        }
    }
    
    Ok(data)
}

/// Extract complete base metadata from SVG circles for metadata-driven architecture
/// CRITICAL: Uses unique ID (node_{position}) for precise extraction, eliminating regex ambiguity
/// This extracts position, base, reactivity, and coordinates (x, y) from SVG
/// CRITICAL FIX: Extracts base type from SVG template's <text> element (actual base) 
/// instead of data-base attribute, ensuring consistency with template
/// The <text> element is found by matching the position in <title> elements
fn extract_base_metadata_from_svg(svg_content: &str) -> Result<Vec<BaseMetadata>, Box<dyn Error>> {
    let mut metadata = Vec::new();
    // Use ID-based extraction for atomic precision (no context window ambiguity)
    let id_regex = Regex::new(r#"id="node_(\d+)""#)?;
    let data_base_regex = Regex::new(r#"data-base="([ATCGU])""#)?; // Fallback: use data-base if text not found
    let reactivity_regex = Regex::new(r#"data-reactivity="([^"]+)""#)?;
    let cx_regex = Regex::new(r#"cx="([^"]+)""#)?;
    let cy_regex = Regex::new(r#"cy="([^"]+)""#)?;
    
    // Build a map of position -> base from <text> elements (actual base in template)
    // CRITICAL: Match by position number from "position.label in template: XXX.Y" format
    // This is the most reliable method as it directly maps actual position to base
    let mut position_to_base: HashMap<u32, char> = HashMap::new();
    let title_template_pos_regex = Regex::new(r"position\.label in template: (\d+)\.([ATCGU])")?; // Actual position and base
    
    // First pass: Extract from "position.label in template: XXX.Y" format (most reliable)
    for line in svg_content.lines() {
        if line.contains("<title>") && !line.contains("numbering-label") {
            if let Some(template_caps) = title_template_pos_regex.captures(line) {
                if let (Ok(actual_pos), Some(base_char)) = (
                    template_caps[1].parse::<u32>(),
                    template_caps[2].chars().next()
                ) {
                    position_to_base.insert(actual_pos, normalize_base(base_char));
                }
            }
        }
    }
    
    // Second pass: Also build coordinate-based map for positions without template info
    // This handles cases where <text> elements exist but don't have "position.label in template" format
    let mut coord_to_base: HashMap<(i32, i32), char> = HashMap::new();
    let title_seq_pos_regex = Regex::new(r"<title>(\d+)")?;
    let text_base_regex = Regex::new(r"<text[^>]*>([ATCGU])</text>")?;
    let text_x_regex = Regex::new(r#"x="([^"]+)""#)?;
    let text_y_regex = Regex::new(r#"y="([^"]+)""#)?;
    
    let lines: Vec<&str> = svg_content.lines().collect();
    for (idx, line) in lines.iter().enumerate() {
        if line.contains("<title>") && !line.contains("numbering-label") {
            // Skip if already processed in first pass
            if title_template_pos_regex.is_match(line) {
                continue;
            }
            
            // Look for <text> element in nearby lines to extract base and coordinates
            if title_seq_pos_regex.is_match(line) {
                for i in (idx.saturating_sub(1))..(idx + 4).min(lines.len()) {
                    if let Some(text_caps) = text_base_regex.captures(lines[i]) {
                        if let Some(base_char) = text_caps[1].chars().next() {
                            // Extract coordinates from the text line
                            if let (Some(x_caps), Some(y_caps)) = (
                                text_x_regex.captures(lines[i]),
                                text_y_regex.captures(lines[i])
                            ) {
                                if let (Ok(x), Ok(y)) = (
                                    x_caps[1].parse::<f64>(),
                                    y_caps[1].parse::<f64>()
                                ) {
                                    // Convert to integer coordinates (multiply by 100 and round) for HashMap key
                                    let int_x = (x * 100.0).round() as i32;
                                    let int_y = (y * 100.0).round() as i32;
                                    coord_to_base.insert((int_x, int_y), normalize_base(base_char));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // CRITICAL FIX: Parse each <circle> tag individually to avoid context window ambiguity
    // Match complete <circle> tags that contain id="node_XXX"
    // Pattern matches both self-closing (<circle ... />) and with closing tag (<circle ...></circle>)
    let circle_with_id_regex = Regex::new(r#"<circle([^>]*?id="node_(\d+)"[^>]*?)(?:/>|>)"#)?;
    
    // Find all circles with unique node ID (physical identifier)
    for circle_cap in circle_with_id_regex.captures_iter(svg_content) {
        // Extract ID from group 2 (the captured number)
        if let Some(id_cap) = circle_cap.get(2) {
            if let Ok(id) = id_cap.as_str().parse::<u32>() {
                // Get the circle tag attributes (group 1 contains all attributes)
                let circle_tag = circle_cap.get(1).unwrap().as_str();
                
                // Extract attributes from this specific circle tag only (no context window)
                let x = cx_regex.captures(circle_tag)
                    .and_then(|c| c[1].parse::<f64>().ok())
                    .unwrap_or(0.0);
                let y = cy_regex.captures(circle_tag)
                    .and_then(|c| c[1].parse::<f64>().ok())
                    .unwrap_or(0.0);
                
                let reactivity = reactivity_regex.captures(circle_tag)
                    .and_then(|c| c[1].parse::<f64>().ok())
                    .unwrap_or(0.0);
                
                // CRITICAL FIX: Extract base from SVG's data-base attribute (most reliable)
                // Parse directly from this circle tag, not from context window
                let base = {
                    // First try: use data-base attribute (set during SVG generation, should be correct)
                    data_base_regex.captures(circle_tag)
                        .and_then(|c| c[1].chars().next())
                        .map(|b| normalize_base(b))
                        .or_else(|| {
                            // Second try: position-based matching (for verification/fallback)
                            position_to_base.get(&id)
                                .copied()
                                .or_else(|| {
                                    // Third try: coordinate-based matching (last resort)
                                    let int_x = (x * 100.0).round() as i32;
                                    let int_y = (y * 100.0).round() as i32;
                                    coord_to_base.get(&(int_x, int_y))
                                        .copied()
                                })
                        })
                }.unwrap_or('N');
                
                // Only add if base is valid
                if base != 'N' {
                    metadata.push(BaseMetadata {
                        id,
                        base,
                        x,
                        y,
                        reactivity,
                    });
                }
            }
        }
    }
    
    Ok(metadata)
}
