# Data Visualization Methods in plot.rs

## Abstract

This document provides a comprehensive description of the data visualization methods implemented in `plot.rs`, including static chart plotting, SVG-based RNA structure visualization, and interactive HTML visualizations. The implementation supports multiple signal types (stop and mutation signals), reactivity data visualization, and advanced interactive features for RNA structure analysis.

## 1. Static Chart Visualization

### 1.1 Signal Distribution Plots

The module provides static chart generation using the `plotters` library for visualizing signal distributions across gene regions.

#### 1.1.1 Stop Signal Distribution Plot

**Function**: `plot_stop_signals()`

**Purpose**: Visualizes the distribution of stop signals (PIPC/Depth) for modified and unmodified samples.

**Input Parameters**:
- `gene_id`: Gene identifier string
- `mod_signals`: Vector of `SignalData` for modified sample
- `unmod_signals`: Vector of `SignalData` for unmodified sample
- `output_dir`: Output directory path

**Output**: PNG image file (`{gene_id}_stop_signals.png`)

**Visualization Details**:
- **Chart Type**: Line plot with scatter points
- **X-axis**: Position (nucleotide position)
- **Y-axis**: Stop Signal (PIPC/Depth)
- **Data Series**:
  - Modified sample: Green line with filled circles (opacity: 0.8, stroke width: 2)
  - Unmodified sample: Blue line with filled circles (opacity: 0.8, stroke width: 2)
- **Chart Dimensions**: 1200 × 800 pixels
- **Legend**: Automatic legend generation with series labels

**Mathematical Formula**:
The stop signal is calculated as:
```
Stop Signal = PIPC / Depth
```
where PIPC (Position-specific Incomplete Primer Count) represents the number of incomplete primers at each position, and Depth represents the sequencing depth.

#### 1.1.2 Mutation Signal Distribution Plot

**Function**: `plot_mutation_signals()`

**Purpose**: Visualizes the distribution of mutation signals for modified and unmodified samples.

**Input Parameters**:
- `gene_id`: Gene identifier string
- `mod_signals`: Vector of `SignalData` for modified sample
- `unmod_signals`: Vector of `SignalData` for unmodified sample
- `output_dir`: Output directory path

**Output**: PNG image file (`{gene_id}_mutation_signals.png`)

**Visualization Details**:
- **Chart Type**: Line plot with scatter points
- **X-axis**: Position (nucleotide position)
- **Y-axis**: Mutation Signal
- **Data Series**:
  - Modified sample: Green line with filled circles
  - Unmodified sample: Blue line with filled circles
- **Chart Dimensions**: 1200 × 800 pixels

**Mathematical Formula**:
The mutation signal represents the mutation rate at each position:
```
Mutation Signal = Mutation Count / Total Depth
```

#### 1.1.3 Reactivity Stop Signal Plot

**Function**: `plot_reactivity_stop_signals()`

**Purpose**: Visualizes reactivity stop signals as color-coded bar charts.

**Input Parameters**:
- `gene_id`: Gene identifier string
- `reactivity_data`: Vector of `ReactivityData`
- `output_dir`: Output directory path

**Output**: PNG image file (`{gene_id}_reactivity_stop.png`)

**Visualization Details**:
- **Chart Type**: Bar chart with color coding
- **X-axis**: Position
- **Y-axis**: Reactivity Stop Signal (0.0 to 1.0)
- **Color Scheme**:
  - Gray (RGB: 128, 128, 128): 0.0 ≤ reactivity ≤ 0.2
  - Yellow (RGB: 255, 255, 0): 0.2 < reactivity ≤ 0.4
  - Red (RGB: 255, 0, 0): 0.4 < reactivity ≤ 1.0
- **Bar Width**: 0.4 units (centered at each position)
- **Legend**: Three color-coded ranges with labels

**Mathematical Formula**:
Reactivity is calculated as:
```
Reactivity = (Modified Signal - Unmodified Signal) / Normalization Factor
```
The reactivity values are normalized to the range [0, 1] for visualization.

#### 1.1.4 Reactivity Mutation Signal Plot

**Function**: `plot_reactivity_mutation_signals()`

**Purpose**: Visualizes reactivity mutation signals as color-coded bar charts.

**Input Parameters**:
- `gene_id`: Gene identifier string
- `reactivity_data`: Vector of `ReactivityData`
- `output_dir`: Output directory path

**Output**: PNG image file (`{gene_id}_reactivity_mutation.png`)

**Visualization Details**:
- **Chart Type**: Bar chart with color coding
- **Color Scheme**: Same as stop signal plot (Gray/Yellow/Red)
- **Bar Width**: 0.4 units

### 1.2 Overall Scatter Plot

**Function**: `plot_overall_scatter()`

**Purpose**: Provides an overview of gene coverage and depth statistics across all analyzed genes.

**Input Parameters**:
- `stats`: HashMap of gene statistics (`GeneStat`)
- `coverage_threshold`: Coverage threshold value (default: 0.2)
- `depth_threshold`: Depth threshold value (default: 50.0)
- `output_path`: Output file path

**Output**: PNG image file (`overall_scatter.png`)

**Visualization Details**:
- **Chart Type**: Scatter plot
- **X-axis**: Coverage (0.0 to max_coverage)
- **Y-axis**: log₁₀(Average Depth + 1)
- **Point Colors**:
  - Red: High coverage and depth (coverage ≥ threshold AND depth ≥ threshold)
  - Gray (RGB: 180, 180, 180): Low coverage or depth
- **Point Size**: 7 pixels
- **Threshold Lines**:
  - Vertical red line: Coverage threshold
  - Horizontal blue line: Depth threshold
- **Chart Dimensions**: 1200 × 900 pixels
- **Legend**: Includes count of high/low groups

**Mathematical Formula**:
```
Y = log₁₀(Average Depth + 1)
```
The logarithmic transformation is applied to handle the wide range of depth values.

## 2. SVG-Based RNA Structure Visualization

### 2.1 Overview

The SVG visualization system maps reactivity data onto RNA secondary structure diagrams. The system supports position alignment, base-specific color mapping, and customizable visualization parameters.

### 2.2 Core Functions

#### 2.2.1 Basic SVG Plotting

**Function**: `plot_reactivity_to_svg()`

**Purpose**: Maps reactivity data to an SVG RNA structure template without reference sequence alignment.

**Input Parameters**:
- `reactivity_file`: CSV file containing reactivity data
- `svg_template`: SVG template file path
- `output_file`: Output SVG file path
- `bases_filter`: String of bases to include (e.g., "ATCG")
- `signal_type`: Signal type ("stop", "mutation", or "all")
- `strand_filter`: Strand filter ("+", "-", or "both")
- `color_ranges`: Optional `BaseColorRanges` structure

**Output**: SVG file with colored circles overlaid on RNA structure

**Data Format**:
The reactivity CSV file should contain columns:
- `position`: Nucleotide position
- `base`: Base type (A, T, C, G, U)
- `reactivity_stop` or `reactivity_mutation`: Reactivity values

#### 2.2.2 SVG Plotting with Alignment

**Function**: `plot_reactivity_to_svg_with_alignment()`

**Purpose**: Maps reactivity data to SVG template with reference sequence alignment to correct position mismatches.

**Input Parameters**:
- `reactivity_file`: CSV file path
- `svg_template`: SVG template file path
- `output_file`: Output SVG file path
- `bases_filter`: Base filter string
- `signal_type`: Signal type
- `strand_filter`: Strand filter
- `ref_sequence_file`: Reference sequence file path (FASTA format)
- `max_shift`: Maximum allowed shift for alignment (default: 5)
- `color_ranges`: Optional color ranges

**Output**: Aligned SVG file

**Alignment Algorithm**:
The alignment process uses a shift-based approach:

1. **Base-specific shift calculation**:
   ```
   shift_base = argmax_shift Σ(position) match_count(CSV_position + shift, SVG_position)
   ```
   where `match_count` counts matching bases between CSV and SVG positions.

2. **Position mapping**:
   ```
   CSV_position = SVG_position + shift_base
   ```

3. **Best shift selection**:
   For each base type, the shift that maximizes the number of matching positions is selected, constrained by `max_shift`.

#### 2.2.3 Color Mapping System

**Structure**: `BaseColorRanges`

The color mapping system uses base-specific reactivity ranges:

```rust
pub struct BaseColorRanges {
    pub a: Vec<ColorRange>,  // Adenine ranges
    pub t: Vec<ColorRange>,  // Thymine/Uracil ranges
    pub c: Vec<ColorRange>,  // Cytosine ranges
    pub g: Vec<ColorRange>,  // Guanine ranges
}

pub struct ColorRange {
    pub min: f64,
    pub max: f64,
    pub color: String,  // RGBA color string
}
```

**Default Color Ranges**:

| Base | Range 1 | Range 2 | Range 3 |
|------|---------|---------|---------|
| A    | 0.0-0.3 (Gray) | 0.3-0.65 (Orange) | 0.65-1.0 (Red) |
| T    | 0.0-0.2 (Gray) | 0.2-0.55 (Orange) | 0.55-1.0 (Red) |
| C    | 0.0-0.3 (Gray) | 0.3-0.6 (Orange) | 0.6-1.0 (Red) |
| G    | 0.0-0.2 (Gray) | 0.2-0.5 (Orange) | 0.5-1.0 (Red) |

**Color Determination Function**: `determine_color()`

**Algorithm**:
```
For each range in color_ranges:
    if is_last_range:
        if reactivity >= range.min AND reactivity <= range.max:
            return range.color
    else:
        if reactivity >= range.min AND reactivity < range.max:
            return range.color

If reactivity < first_range.min:
    return first_range.color
If reactivity > last_range.max:
    return last_range.color
Else:
    return gray (no data)
```

**Mathematical Formula**:
The color assignment follows a piecewise function:
```
C(reactivity, base) = {
    color₁  if reactivity ∈ [0, threshold₁)
    color₂  if reactivity ∈ [threshold₁, threshold₂)
    color₃  if reactivity ∈ [threshold₂, 1.0]
}
```
where `threshold₁` and `threshold₂` are base-specific values.

#### 2.2.4 SVG Drawing Configuration

**Structure**: `SvgDrawConfig`

**Parameters**:
- `circle_radius`: Radius of reactivity circles (default: 2.8)
- `circle_stroke_width`: Stroke width for circles (default: 0.2)
- `circle_has_stroke`: Whether to draw circle stroke (default: false)
- `circle_filled`: Whether circles are filled or hollow (default: false)
- `font_color`: Optional font color override
- `circles_before_text`: Layer ordering flag (default: false, text on top)
- `draw_color_bars`: Whether to draw color legend bars (default: true)
- `legend_item_width`: Width of legend items for horizontal layout (default: 30.0)
- `legend_item_height`: Height of legend items for vertical layout (default: 15.0)

**Circle Rendering**:
For filled circles:
```
style = "fill: {color}; stroke: {color}; stroke-width: {stroke_width};"
```

For hollow circles:
```
style = "fill: none; stroke: {color}; stroke-width: {stroke_width};"
```

## 3. Interactive HTML Visualization

### 3.1 Overview

The interactive HTML visualization system provides a comprehensive web-based interface for exploring RNA structure data with reactivity overlays. It includes advanced features such as zoom, pan, filtering, tooltips, and dynamic color adjustment.

**Function**: `generate_interactive_html()`

**Input Parameters**:
- `svg_content`: SVG content string (with circles already drawn)
- `reactivity_data`: Vector of `SvgReactivityData`
- `output_file`: Output HTML file path
- `signal_type`: Signal type identifier
- `color_ranges`: Initial `BaseColorRanges` configuration
- `initial_circle_filled`: Initial circle fill state

**Output**: Self-contained HTML file with embedded SVG and JavaScript

### 3.2 Core Interactive Features

#### 3.2.1 Zoom and Pan

**Zoom Implementation**:
- **Mouse Wheel Zoom**: Zoom in/out centered on mouse cursor position
- **Zoom Formula**:
  ```
  delta = (deltaY > 0) ? 1.1 : 0.9
  mouse_x = (clientX - rect.left) / rect.width * viewBox.width + viewBox.x
  mouse_y = (clientY - rect.top) / rect.height * viewBox.height + viewBox.y
  viewBox.x = mouse_x - (mouse_x - viewBox.x) * delta
  viewBox.y = mouse_y - (mouse_y - viewBox.y) * delta
  viewBox.width *= delta
  viewBox.height *= delta
  ```

- **Keyboard Shortcuts**:
  - `+` or `=`: Zoom in
  - `-` or `_`: Zoom out
  - `0` or `F`: Fit to view
  - `R`: Reset view

**Pan Implementation**:
- **Methods**:
  - Middle mouse button drag
  - Spacebar + left mouse button drag
  - Hand tool mode (toolbar button)
- **Pan Formula**:
  ```
  scaleX = viewBox.width / rect.width
  scaleY = viewBox.height / rect.height
  deltaX = (panStart.x - clientX) * scaleX
  deltaY = (panStart.y - clientY) * scaleY
  viewBox.x = panStart.viewBoxX + deltaX
  viewBox.y = panStart.viewBoxY + deltaY
  ```

#### 3.2.2 Minimap Navigation

**Purpose**: Provides an overview of the entire structure with a viewport indicator.

**Implementation**:
- **Visibility**: Shown when zoomed in (viewBox.width < 0.9 × initialViewBox.width)
- **Viewport Calculation**:
  ```
  scaleX = 100 / initialViewBox.width
  scaleY = 100 / initialViewBox.height
  viewport.x = (viewBox.x - initialViewBox.x) * scaleX
  viewport.y = (viewBox.y - initialViewBox.y) * scaleY
  viewport.width = viewBox.width * scaleX
  viewport.height = viewBox.height * scaleY
  ```
- **Click Navigation**: Clicking on minimap centers viewport at clicked position

#### 3.2.3 Tooltip System

**Function**: `showTooltip()`

**Display Information**:
- Position (nucleotide position)
- Base type (A, T, C, G, U)
- Reactivity value (formatted to 4 decimal places)
- Signal type

**Position Calculation**:
```
tooltip.left = clientX + 10px
tooltip.top = clientY + 10px
```

**Trigger**: Mouse hover over reactivity circles

#### 3.2.4 Filtering System

**Reactivity Threshold Filtering**:
- **Min/Max Sliders**: Filter circles by reactivity range
- **Implementation**:
  ```
  if circle.reactivity < minThreshold OR circle.reactivity > maxThreshold:
      circle.style.display = 'none'
  else:
      circle.style.display = ''
  ```

**Base Type Filtering**:
- **Checkboxes**: Enable/disable display of specific base types (A, T, C, G)
- **Implementation**:
  ```
  if baseEnabled[base] == false:
      circle.style.display = 'none'
  ```

#### 3.2.5 Dynamic Color Range Adjustment

**Unified Mode**:
- Single set of thresholds and colors applied to all bases
- **Controls**:
  - `unifiedRange1`: Low-Mid threshold (0.0-1.0)
  - `unifiedRange2`: Mid-High threshold (0.0-1.0)
  - `unifiedColor1`: Low reactivity color
  - `unifiedColor2`: Mid reactivity color
  - `unifiedColor3`: High reactivity color

**Individual Mode**:
- Separate thresholds and colors for each base type
- **Controls per base**:
  - `range{Base}1`: Low-Mid threshold
  - `range{Base}2`: Mid-High threshold
  - `color{Base}1`: Low color
  - `color{Base}2`: Mid color
  - `color{Base}3`: High color

**Color Update Formula**:
```
For each base:
    colorRanges[base] = [
        {min: 0.0, max: threshold1, color: color1},
        {min: threshold1, max: threshold2, color: color2},
        {min: threshold2, max: 1.0, color: color3}
    ]
```

**Real-time Application**:
Colors are recalculated and applied immediately when thresholds or colors change using the `applyColorAndVisibility()` function.

#### 3.2.6 Selection and Editing Tools

**Tool Modes**:
1. **Select Mode**: Click circles to show context menu for individual editing
2. **Box Select Mode**: Drag to select multiple bases for batch editing
3. **Hand Mode**: Pan tool for navigation

**Context Menu (Select Mode)**:
- Edit font size
- Edit font color
- Edit circle radius
- Edit circle color
- Reset base style

**Batch Edit Menu (Box Select Mode)**:
- Batch edit font size
- Batch edit font color
- Batch edit circle radius
- Batch edit circle color
- Batch reset style

**Selection Box Calculation**:
```
selectedBases = {position | circle.getBoundingClientRect() intersects selectionBox}
```

#### 3.2.7 Legend Customization

**Legend Parameters**:
- **Position**: Top, Bottom, Left, Right
- **Direction**: Horizontal or Vertical
- **Font Size**: Adjustable (default: 12px)
- **Item Dimensions**:
  - Width: For horizontal layout (default: 30px)
  - Height: For vertical layout (default: 15px)

**Legend Layout Calculation**:

For horizontal direction:
```
legendWidth = itemsCount × itemWidth
legendHeight = itemHeight + textHeight + spacing
```

For vertical direction:
```
legendWidth = max(itemWidth, textWidth)
legendHeight = itemsCount × itemHeight + textHeight + spacing
```

**Position Calculation** (example for bottom position):
```
baseX = (svgWidth / 2) - (legendWidth / 2)
baseY = svgHeight - totalLegendsHeight - 20
```

#### 3.2.8 Export Functionality

**Supported Formats**:
1. **SVG Export**:
   - Direct serialization of SVG element
   - Preserves all styling and structure
   - Formula: `serializer.serializeToString(svg)`

2. **PNG Export**:
   - Rasterization of SVG to canvas
   - **Resolution Control**: Multiplier (default: 1.0)
   - **DPI Control**: Dots per inch (default: 72)
   - **Canvas Dimensions**:
     ```
     canvas.width = img.width × resolution × (dpi / 72)
     canvas.height = img.height × resolution × (dpi / 72)
     ```
   - Formula: `canvas.toBlob(blob, 'image/png')`

**Export Parameters**:
- File name: User-specified (default: "rna_structure_interactive")
- Format: SVG or PNG
- PNG resolution: 0.5x to 5.0x
- PNG DPI: 72 to 600

#### 3.2.9 Undo/Redo System

**Implementation**:
- **State Storage**: Saves circle and text element states
- **State Structure**:
  ```javascript
  {
      circles: [{position, r, style}],
      texts: [{text, fontSize, fill}]
  }
  ```
- **Stack Management**:
  - Maximum undo steps: 50
  - Redo stack cleared on new action
- **Keyboard Shortcuts**:
  - `Ctrl+Z` or `Cmd+Z`: Undo
  - `Ctrl+Shift+Z` or `Cmd+Shift+Z`: Redo
  - `Ctrl+Y` or `Cmd+Y`: Redo

**State Save Triggers**:
- Font size/color changes
- Circle radius/color changes
- Style resets
- Batch edits

#### 3.2.10 Search and Information Panel

**Search Functionality**:
- **Keyboard Shortcut**: `S`
- **Features**:
  - Search by position
  - Search by base type
  - Highlight matching positions
  - Navigate to position

**Information Panel**:
- **Keyboard Shortcut**: `I`
- **Display**:
  - Total positions
  - Visible positions
  - Filtered positions
  - Current zoom level
  - Viewport coordinates

**Update Frequency**: Every 2 seconds

### 3.3 Data Structure

**SvgReactivityData**:
```rust
pub struct SvgReactivityData {
    pub position: u32,      // Nucleotide position
    pub base: char,         // Base type (A, T, C, G, U)
    pub reactivity: f64,    // Reactivity value (0.0-1.0)
}
```

**JavaScript Data Format**:
```javascript
reactivityData = {
    signalType: "stop" | "mutation" | "all",
    positions: {
        [position]: {
            base: "A" | "T" | "C" | "G" | "U",
            reactivity: 0.0-1.0
        }
    }
}
```

### 3.4 Performance Optimizations

1. **Event Listener Management**: Clones SVG elements to remove old listeners
2. **Debounced Updates**: 500ms delay for state saving
3. **Conditional Rendering**: Only updates visible elements
4. **ViewBox-based Zoom**: Efficient SVG viewport manipulation

## 4. Multi-Signal Visualization

### 4.1 Automatic Signal Detection

**Function**: `detect_signal_types()`

**Purpose**: Automatically detects available signal types from CSV file headers.

**Algorithm**:
```
1. Read first line of CSV file
2. Split by comma delimiter
3. For each field:
    if field contains "reactivity_":
        signal_name = field.replace("reactivity_", "")
        add signal_name to signal_types
4. If no reactivity columns found:
    add "score" as default signal type
```

**Function**: `plot_multiple_signals_to_svg()`

**Purpose**: Generates separate SVG files for each detected signal type.

**Output Files**:
- `rna_structure_colored_{signal_type}.svg` for each signal type

## 5. Implementation Details

### 5.1 Threading and Parallelization

The `plot_signal_distributions()` function uses Rayon for parallel processing:

```rust
plot_tasks.into_par_iter()
    .map(|task| {
        // Plot generation for each gene
    })
    .collect()
```

**Thread Configuration**:
- Default: Auto-detect CPU cores
- Custom: User-specified thread count

### 5.2 File I/O

**Input Formats**:
- CSV: Comma-separated values with header row
- SVG: Standard SVG format with `<title>` elements for positions
- FASTA: Reference sequence for alignment

**Output Formats**:
- PNG: Bitmap images for static charts
- SVG: Vector graphics for structure visualization
- HTML: Self-contained interactive visualizations

### 5.3 Error Handling

- File not found: Returns error with file path
- Invalid data format: Skips invalid entries with warnings
- Alignment failures: Falls back to non-aligned mode
- Rendering errors: Logs error and continues with next gene

## 6. Usage Examples

### 6.1 Static Chart Generation

```rust
plot_signal_distributions(
    "modified.csv",
    "unmodified.csv",
    "output/",
    Some(8),           // 8 threads
    Some(0.2),          // Coverage threshold
    Some(50.0),        // Depth threshold
    Some("reactivity.csv"),
    Some("annotation.gff"),
    &mut logger
)?;
```

### 6.2 SVG Visualization

```rust
plot_reactivity_to_svg(
    "reactivity.csv",
    "structure.svg",
    "output_colored.svg",
    "ATCG",
    "stop",
    "+",
    None  // Use default color ranges
)?;
```

### 6.3 Interactive Visualization

```rust
generate_interactive_visualization(
    "reactivity.csv",
    "structure_colored.svg",
    "interactive.html",
    "ATCG",
    "stop",
    "+",
    None,  // Use default color ranges
    false  // Hollow circles
)?;
```

## 7. References

- Plotters Library: https://github.com/38/plotters
- SVG Specification: https://www.w3.org/TR/SVG2/
- Spectrum Color Picker: https://seballot.github.io/spectrum/

## 8. Conclusion

The visualization system in `plot.rs` provides a comprehensive suite of tools for RNA structure analysis, ranging from static statistical charts to highly interactive web-based visualizations. The implementation supports multiple signal types, flexible color mapping, and advanced user interaction features, making it suitable for both automated analysis pipelines and interactive exploration of RNA modification data.






