# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [0.15.4] - 2026-02-08

### Fixed
- **Base Mismatch in SVG Metadata Extraction**: Fixed critical bug where `extract_reactivity_from_svg()` and `extract_base_metadata_from_svg()` functions were matching attributes from adjacent circles due to context window method
  - **Issue**: SVG circles are often on the same line or adjacent blocks, and the context window method (200 characters around ID match) was capturing attributes from previous circles, causing base type mismatches
  - **Example**: `node_498` (data-base="C") was getting base="G" from adjacent `node_584` (data-base="G")
  - **Fix**: Replaced context window method with individual circle tag parsing - each `<circle>` tag is now parsed independently to extract only its own attributes
  - **Impact**: Eliminates "Base mismatch" errors in HTML metadata verification, ensures `RNA_METADATA` matches circle `data-base` attributes correctly

- **Base Type Extraction from SVG Template**: Fixed base type extraction in `parse_svg_template_for_circles()` to use reliable template position mapping
  - **Issue**: Previous code used "nearby text search" method which could match wrong `<text>` elements
  - **Fix**: Added `extract_template_pos_to_base()` function to extract base mapping from `position.label in template: POS.BASE` format (most reliable source)
  - **Implementation**: `parse_svg_template_for_circles()` now uses `template_pos2base` mapping instead of searching nearby lines
  - **Impact**: Ensures `data-base` attribute uses SVG template's actual base type (`normalized_svg_base`) instead of CSV base type

### Changed
- **SVG Circle Parsing**: Modified `extract_reactivity_from_svg()` to parse each `<circle>` tag individually
  - Uses regex to match complete circle tags (self-closing or with closing tag)
  - Extracts attributes (`data-position`, `data-base`, `data-reactivity`) only from the matched circle tag
  - Eliminates cross-contamination between adjacent circles

- **SVG Metadata Extraction**: Modified `extract_base_metadata_from_svg()` to parse each `<circle>` tag individually
  - Uses regex to match circles with `id="node_XXX"` pattern
  - Extracts attributes (`cx`, `cy`, `data-base`, `data-reactivity`) only from the matched circle tag
  - Maintains fallback to position-based and coordinate-based matching for verification

### Technical Details
- **File**: `src/plot.rs`
  - `extract_reactivity_from_svg()`: Changed from context window to individual circle tag parsing using regex pattern `<circle[^>]*?/>|<circle([^>]*?)>(?:[^<]*</circle>)?`
  - `extract_base_metadata_from_svg()`: Changed from context window to individual circle tag parsing using regex pattern `<circle([^>]*?id="node_(\d+)"[^>]*?)(?:/>|>)`
  - `extract_template_pos_to_base()`: New function to extract `position -> base` mapping from `position.label in template: POS.BASE` format
  - `parse_svg_template_for_circles()`: Updated to accept `template_pos2base` parameter and use it instead of nearby text search
  - `draw_svg_with_custom_layers()`: Updated to call `extract_template_pos_to_base()` and pass mapping to `parse_svg_template_for_circles()`

### Impact
- **Metadata Accuracy**: `RNA_METADATA` in generated HTML now correctly matches circle `data-base` attributes
- **Verification**: `Metadata verification` should now pass without "Base mismatch" errors
- **Threshold Control**: Eliminates cross-contamination where adjusting one base's threshold affects other bases
- **Backward Compatible**: Fully backward compatible; no API changes, only internal parsing improvements

# [0.14.10] - 2026-02-03

### Changed
- **Zarringhalam remap negative value handling**: Integrated the "align min to 0" logic into `zarringhalam_remap` internal implementation
  - **Previously**: In the norm workflow, rescale was done separately before remap. If only remap was used without norm, negative values would still be lost
  - **Now**: When using Zarringhalam mapping, `zarringhalam_remap` internally first performs linear scaling with the current sequence minimum as 0 ([min, 1] → [0, 1]), then applies piecewise linear mapping
  - **Effect**: Regardless of whether norm is performed first, as long as Zarringhalam is used, it will preserve mod/unmod relative differences, avoiding the shortcomings of Siegfried's negative value information loss and amplification of unmod signal imbalance
  - **Implementation**: In `src/norm.rs`, `zarringhalam_remap` first calculates `min_val` and `range = 1 - min_val`, then applies `(x - min_val) / range` to the input before applying Zarringhalam piecewise linear mapping; degenerate cases (range ≤ 1e-10) are mapped to 0

### Technical Details
- **File**: `src/norm.rs`
  - `zarringhalam_remap`: Internally first rescales [min, 1] → [0, 1], then executes Zarringhalam piecewise linear mapping
  - Removed independent `rescale_min_to_zero` and separate calls in the norm workflow, only calls `zarringhalam_remap` once when `--linear` is selected

# [0.14.9] - 2026-02-01

### Fixed
- **Interactive SVG Right-Click Menu Style Changes**: Fixed issue where style changes from right-click context menu were not taking effect
  - **Issue**: Right-click menu changes to font size, font color, circle radius, and circle color were not being applied to SVG elements
  - **Fix**: Improved DOM manipulation to use both `setAttribute()` and direct `style` object updates for reliable style application
  - **Implementation**: 
    - `editBaseFontSize()`: Added input validation and ensures both `setAttribute('font-size')` and `text.style.fontSize` are updated
    - `editBaseFontColor()`: Added hex color format validation and ensures both `setAttribute('fill')` and `text.style.fill` are updated
    - `editBaseCircleRadius()`: Added input validation and ensures `setAttribute('r')` is updated
    - `editBaseCircleColor()`: Replaced unreliable regex-based style manipulation with direct DOM `style` object manipulation, correctly detects filled vs hollow circles, preserves `stroke-width`, and updates both `style` and `setAttribute` for consistency
    - `resetBaseStyle()`: Improved to reset font size, font color, and circle radius by updating both `setAttribute` and `style` properties
  - **Impact**: Users can now successfully modify element styles through the right-click context menu

- **Text Element Not Found Error**: Fixed "Text element not found for this position" error when editing font styles
  - **Issue**: `findTextForPosition()` function could not locate text elements for some positions, preventing font style editing
  - **Fix**: Enhanced `findTextForPosition()` with multiple fallback methods:
    - Method 1: Checks for text element within the same `<g>` group as the circle
    - Method 2: Iterates through all `<g>` elements, looking for `<title>` tags that match the position (handles formats like "790" or "790 (position.label in template: ...)")
    - Method 3: Checks sibling `<g>` elements as fallback
  - **Impact**: Font size and color editing now works reliably for all positions

- **State Saving Null Reference Errors**: Fixed `TypeError: Cannot read properties of null (reading 'value')` errors in `saveState()` function
  - **Issue**: `saveState()` attempted to access DOM element properties without checking if elements exist
  - **Fix**: Added null checks for all DOM element accesses before reading `.value` or `.checked` properties
  - **Impact**: Prevents JavaScript errors when saving application state

### Changed
- **Legend Initial Position**: Changed legend initial position to avoid obscuring RNA secondary structure
  - **Previous**: Legend was positioned at top-left corner, often obscuring the RNA structure
  - **New**: Legend is now positioned at bottom-right corner of the SVG canvas
  - **Implementation**: 
    - Added `extract_svg_width()` and `extract_svg_height()` functions to extract SVG canvas dimensions from `viewBox` or `width`/`height` attributes
    - Modified `draw_color_bar()` to calculate legend position based on SVG dimensions: `base_x = (svg_width - legend_width - margin).max(margin)` and `base_y = (svg_height - total_legend_height - margin) + position * row_height`
  - **Impact**: Legend no longer obscures the RNA structure visualization

- **Page Loading Experience**: Improved initial page loading experience with better loading message
  - **Previous**: Brief flash of "drag to upload" message before SVG loads
  - **New**: Shows "Loading xxx.svg..." message for 2.5 seconds before displaying visualization
  - **Implementation**:
    - Extracts SVG filename from output file path (e.g., "EC_16S.svg" from HTML filename)
    - Displays loading spinner with filename
    - Upload area is hidden during loading, showing only the loading message
  - **Impact**: Users see a clear loading indication instead of a brief flash

- **Minimap Viewport Dragging**: Enabled dragging of minimap viewport (focus box) to navigate the main view
  - **Previous**: Viewport was locked in center, could only click minimap background to jump to position
  - **New**: Viewport can be dragged to move the main view in real-time
  - **Implementation**:
    - Added `mousedown`, `mousemove`, and `mouseup` event listeners to viewport element using capture phase
    - Calculates viewport position in minimap coordinates (0-100) and updates main SVG `viewBox` accordingly
    - Constrains viewport movement to minimap bounds
  - **Impact**: Users can now drag the viewport to navigate the visualization more intuitively

- **Minimap Auto-Hide Behavior**: Improved minimap auto-hide timing to only hide after user activity stops
  - **Previous**: Minimap would hide 3 seconds after being displayed, even during active user operations
  - **New**: Minimap hides 3 seconds after the last user activity (zoom, pan, or viewport drag)
  - **Implementation**:
    - Added debounce mechanism: 100ms delay after last activity before starting 3-second hide timer
    - Timer is reset on any user activity: zoom operations (`zoomIn()`, `zoomOut()`, mouse wheel), pan operations, viewport dragging, or minimap background clicks
    - `resetMinimapHideTimer()` function manages the debounce and hide timer logic
  - **Impact**: Minimap stays visible during active use and only hides when user is idle

- **Version Number Display**: Changed HTML version number to automatically use version from `Cargo.toml`
  - **Previous**: Version number was hardcoded in HTML generation code (`v0.11.6`)
  - **New**: Version number is read from `crate::VERSION` constant (defined in `main.rs` using `env!("CARGO_PKG_VERSION")`)
  - **Implementation**: 
    - Modified `generate_interactive_html()` to use `crate::VERSION` instead of hardcoded string
    - Version is inserted into HTML header and footer using `html.push_str(version)`
  - **Impact**: Version number in HTML output automatically matches the project version in `Cargo.toml`, no manual updates needed

### Technical Details
- **File**: `src/plot.rs`
  - Enhanced JavaScript functions for right-click menu style editing with proper DOM manipulation
  - Improved `findTextForPosition()` with multiple fallback methods for robust text element location
  - Added `extract_svg_width()` and `extract_svg_height()` functions for SVG dimension extraction
  - Modified `draw_color_bar()` to position legend at bottom-right based on SVG dimensions
  - Updated page loading logic to show loading message with extracted SVG filename
  - Implemented minimap viewport dragging with event capture phase for reliable event handling
  - Added minimap auto-hide timer with debounce mechanism for better user experience
  - Changed version number source from hardcoded string to `crate::VERSION` constant

# [0.14.8] - 2026-01-30

### Added
- **Bedgraph Format Support in `modtector convert`**: Implemented bedgraph file format conversion directly in Rust
  - New `--format bedgraph` option for `modtector convert` command
  - Supports bedgraph format: `chr\tstart\tend\tcoverage\tstrand`
  - For 5' end counting (bedtools genomecov -5), coverage represents depth (total reads at that position)
  - RT stop count is calculated as `coverage[i] - coverage[i+1]` for consecutive positions on the same strand
  - Automatic format detection: Detects bedgraph format by checking for 5 tab-separated fields with numeric start/end/coverage and strand (+/-)
  - Optional `--ref-fasta` parameter: If provided, extracts reference base information; otherwise sets base to 'N'
  - Supports `--filter-strand` parameter to filter by strand orientation
  - Replaces Python script-based conversion in workflow with native Rust implementation for better performance and Python version compatibility

### Fixed
- **Bedgraph Depth Calculation**: Fixed incorrect depth calculation in bedgraph format conversion
  - **Issue**: Previously treated bedgraph `coverage` as RT stop count, causing `depth = stop_count` and resulting in 100% stop signal ratios
  - **Fix**: Correctly interprets `coverage` as depth (total reads), and calculates RT stop count from adjacent position differences
  - **Impact**: Stop signal ratios are now correctly calculated (e.g., ~6.6% positions with stop/depth > 0.9 instead of 100%)
  - **Implementation**: Two-pass algorithm that first reads all bedgraph entries, groups by (chr, strand), sorts by position, then calculates RT stop counts using consecutive position coverage differences

### Changed
- **Workflow Integration**: Updated `runs/08-MaP/Snakefile`
  - `count_bt` rule now correctly finds and renames bedgraph output files from `call_stop.py`
  - `convert_bt_to_pileup` rule now uses `modtector convert` instead of embedded Python script
  - Improved performance and maintainability by using native Rust implementation
  - Better error handling and logging through modtector's built-in logging system

### Technical Details
- **File**: `src/convert.rs`
  - Implemented `convert_bedgraph_to_pileup` function for bedgraph file parsing and conversion
  - Two-pass algorithm: First pass reads and groups entries by (chr, strand), second pass calculates RT stop counts
  - Updated `detect_input_format` to recognize bedgraph format (5 tab-separated fields: chr, start, end, coverage, strand)
  - Converts 0-based start position to 1-based position for pileup format
  - Handles strand filtering, reference base extraction (if reference provided), and error reporting
  - Integrated with existing `load_all_reference_sequences` function for multi-sequence reference support

### Usage
```bash
# Convert bedgraph file to pileup format (with reference for base information)
modtector convert \
    -i input.bedgraph \
    -o output.csv \
    -f bedgraph \
    --ref-fasta reference.fa \
    --filter-strand "+"

# Convert without reference (bases will be 'N')
modtector convert \
    -i input.bedgraph \
    -o output.csv \
    -f bedgraph \
    --filter-strand "+"

# Process all strands
modtector convert \
    -i input.bedgraph \
    -o output.csv \
    -f bedgraph \
    --ref-fasta reference.fa \
    --filter-strand None
```

# [0.14.7] - 2026-01-30

### Added
- **Native icSHAPE RT Format Support in `modtector convert`**: Implemented icSHAPE-pipe RT file format conversion directly in Rust
  - New `--format icSHAPE-rt` option for `modtector convert` command
  - New `--filter-strand` parameter: Filter by strand orientation ('+' or '-' or 'None' to process all strands, default: '+')
  - Automatic format detection: Detects icSHAPE RT format by checking for `@ColNum` or `@ChrID` header lines
  - Requires `--ref-fasta` parameter to load reference sequences for base information extraction
  - Replaces Python script-based conversion in workflow with native Rust implementation for better performance

### Changed
- **Workflow Integration**: Updated `runs/08-MaP/Snakefile`
  - `convert_ic_to_pileup` rule now uses `modtector convert` instead of embedded Python script
  - Improved performance and maintainability by using native Rust implementation
  - Better error handling and logging through modtector's built-in logging system

### Technical Details
- **File**: `src/convert.rs`
  - Added `filter_strand: Option<String>` field to `ConvertArgs` struct
  - Implemented `convert_icshape_rt_to_pileup` function for RT file parsing and conversion
  - Updated `detect_input_format` to recognize icSHAPE RT format (starts with `@ColNum` or `@ChrID`)
  - Integrated with existing `load_all_reference_sequences` function for multi-sequence reference support
  - Handles comment lines (starting with `@`), strand filtering, and reference base extraction

### Usage
```bash
# Convert icSHAPE RT file to pileup format
modtector convert \
    -i input.rt \
    -o output.csv \
    -f icSHAPE-rt \
    --ref-fasta reference.fa \
    --filter-strand "+"

# Process all strands
modtector convert \
    -i input.rt \
    -o output.csv \
    -f icSHAPE-rt \
    --ref-fasta reference.fa \
    --filter-strand None
```

# [0.14.6] - 2026-01-30

### Added
- **icSHAPE RT Format Support for `modtector convert`**: Extended `convert` command to support icSHAPE-pipe RT file format conversion
  - Support for converting icSHAPE-pipe RT files (tab-separated format: `chr_id\tstrand\tposition\tRT_count\tBD_count`) to modtector pileup CSV format
  - New `filter_strand` parameter in workflow: Allows filtering by strand orientation (e.g., only process '+' strand) when converting icSHAPE RT files
  - Automatically skips comment lines (starting with `@`) in RT files
  - Handles both positive and negative strand data based on parameter configuration
  - Enables integration of icSHAPE-pipe workflow results into modtector analysis pipeline

### Technical Details
- **Workflow Integration**: `runs/08-MaP/Snakefile`
  - Added `convert_ic_to_pileup` rule for converting icSHAPE RT format to pileup
  - RT file format: `chr_id\tstrand\tposition\tRT_count\tBD_count`
  - Pileup output includes: RT stops count (`pipe_truncation_count`), depth (RT_count + BD_count), and reference base information
  - Configurable strand filtering via `filter_strand` parameter (default: "+" for single-direction reference)

### Usage
```bash
# In Snakemake workflow, convert icSHAPE RT file to pileup format
# filter_strand parameter controls which strand to process:
# - "+" or "-": only process specified strand
# - None: process all strands
```

# [0.14.5] - 2026-01-30

### Added
- **Dual Input Mode for `modtector convert`**: Extended `convert` command to support simultaneous input of mutation and stop signal files
  - New `--input-mutation` parameter: Optional mutation input file (rf-rctools format generated with `rf-count -m`)
  - New `--input-stop` parameter: Optional stop input file (rf-rctools format generated without `-m` flag)
  - When both files are provided, they are merged into a single pileup output with both mutation and stop counts
  - Missing files are handled gracefully: corresponding columns (mutation_count or pipe_truncation_count) are set to 0
  - Enables processing of datasets that have both mutation and stop signals from RNA Framework

### Technical Details
- **File**: `src/convert.rs`
  - Added `input_mutation` and `input_stop` optional fields to `ConvertArgs` struct
  - Created `merge_rf_rctools_to_pileup` function to merge two rf-rctools files
  - Updated `validate_convert_args` to allow missing files in dual input mode
  - Modified `convert_bamreadcount_to_pileup` to detect and handle dual input mode
  - Merge function combines data by (chr, position) keys, preserving both mutation and stop counts

### Usage
```bash
# Merge mutation and stop files
modtector convert --input-mutation mutation.tsv --input-stop stop.tsv -o output.csv -f rf-rctools

# Use only mutation file (stop columns will be 0)
modtector convert --input-mutation mutation.tsv -o output.csv -f rf-rctools

# Use only stop file (mutation columns will be 0)
modtector convert --input-stop stop.tsv -o output.csv -f rf-rctools
```

# [0.14.4] - 2025-01-XX

### Added
- **Quality-Based Effective Depth Calculation**: Enhanced depth filtering based on base quality scores
  - New `effective_depth` field in output CSV files when `--min-base-qual` is enabled
  - `effective_depth` counts only bases with quality >= threshold (similar to shapemapper2)
  - `depth` field still reports total read depth (all bases) for comparison
  - Quality filtering improves mutation rate calculation accuracy by excluding low-quality bases

- **PCR Bias Correction Method**: Chi-Square distribution-based PCR amplification bias correction
  - Implemented in `modtector correct` command (available but disabled in default workflow)
  - Uses Chi-Square distribution fitting to model depth-mutation rate relationship
  - Applies correction factors to adjust effective depth based on PCR bias
  - Configurable weights: `--weight-increase` and `--weight-decrease` parameters
  - Can be combined with quality filtering for comprehensive depth correction

### Changed
- **Snakefile Workflow**: Updated to use quality-filtered depth only
  - `rule count_pileup`: Uses `--min-base-qual 20` to generate `effective_depth` column
  - `rule correct_pileup`: Disabled (commented out) - PCR bias correction optional
  - `rule reactivity`: Now uses quality-filtered `pileup.csv` files directly
  - All downstream rules updated to use `pileup.csv` instead of `pileup_corrected.csv`

- **Output Format**: CSV files now include `effective_depth` column when quality filtering is enabled
  - Format: `ChrID,Strand,ChrPos,Base,Count,StopCount,depth,effective_depth,ins,del,A,C,G,T`
  - `effective_depth` column appears between `depth` and `ins` columns
  - Backward compatible: `effective_depth` column only added when `--min-base-qual > 0` or `--pcr-bias-correction` enabled

### Technical Details
- **File**: `src/main.rs`
  - Added `effective_depth: usize` field to `Item` struct
  - Quality filtering logic in `process_region_chunk_ultra_fast` function
  - Output logic updated to include `effective_depth` column conditionally
- **File**: `src/correct.rs`
  - PCR bias correction implementation (Chi-Square fitting)
  - Correction factor calculation based on depth-mutation rate relationship
- **File**: `runs/08-MaP/Snakefile`
  - Updated workflow to use quality filtering by default
  - PCR bias correction rule commented out but preserved

### Impact
- **Improved Accuracy**: Quality-filtered depth provides more accurate mutation rate calculation
  - Excludes low-quality bases that may introduce noise
  - Similar approach to shapemapper2's `effective_depth` calculation
- **Workflow Simplification**: Using quality filtering alone simplifies the pipeline
  - Single-step depth correction (quality filtering) instead of two-step (quality + PCR correction)
  - PCR bias correction available as optional advanced feature
- **Backward Compatible**: Existing workflows continue to work
  - Quality filtering disabled by default (`--min-base-qual 0`)
  - `effective_depth` column only added when filtering enabled

### References
- Detailed documentation: `note/v0.14.4-work-DEPTH_CORRECTION_METHODS.md`
- shapemapper2 reference: `MutationProcessing.h` - `filterQscoresCountDepths` function
- PCR bias correction: Based on Chi-Square distribution model

# [0.14.0] - 2025-01-XX

### Added
- **Base Quality Filtering for Mutation Detection**: New per-base quality filtering in `modtector count` command
  - New `--min-base-qual` parameter (default: 0, disabled) to filter mutations by base quality score
  - Only mutations with quality >= threshold are counted (per-base filtering like RNAFramework rf-count)
  - Recommended value: 20 (same as RNAFramework rf-count)
  - Quality scores are decoded Phred scores (0-93), not ASCII
  - Significantly improves signal-to-noise ratio by filtering low-quality mutations

### Changed
- **Mutation Counting Logic**: Enhanced `modtector count` to support base quality filtering
  - Quality check only applied to mutations (not reference matches)
  - Filtering disabled by default (min-base-qual=0) for backward compatibility
  - Quality filtering applied per-base using Phred scores from BAM file
  - Updated Snakefile `rule count_pileup` to use `--min-base-qual 20` by default

### Fixed
- **Critical Bug Fix**: Fixed incorrect quality score comparison logic
  - Previous code incorrectly added 33 to threshold (treating as Phred+33 ASCII encoding)
  - Now correctly compares decoded Phred scores directly
  - Fixes issue where all mutations were filtered when quality filtering was enabled

### Impact
- **Signal Quality**: Base quality filtering filters ~15% of low-quality mutations
  - Example: 377,992 mutations → 318,891 mutations (15% filtered) with min-base-qual=20
  - Improves signal-to-noise ratio by removing unreliable mutation calls
  - Better alignment with RNAFramework rf-count behavior
- **Performance**: Minimal performance impact (~0.9% overhead)
  - Quality check only performed for mutations, not all bases
  - Efficient per-base quality score lookup from BAM records
- **Backward Compatible**: Fully backward compatible; filtering disabled by default
  - Existing workflows continue to work without modification
  - New filtering only activated with explicit `--min-base-qual` parameter
  - Snakefile updated to enable filtering for new runs

### Technical Details
- **File**: `src/main.rs`
- **Function**: `process_region_chunk_ultra_fast` (lines 3162-3174)
- **Quality Score Handling**: 
  - Uses `rust_htslib` `record.qual()` which returns decoded Phred scores (0-93)
  - Direct comparison: `base_qual >= min_base_qual` (no offset needed)
  - Quality check only for mutations: `is_mutation && min_base_qual > 0`
- **Snakefile Integration**: Updated `rule count_pileup` in `runs/08-MaP/Snakefile`
  - Added `--min-base-qual 20` parameter to modtector count command
  - All md method count outputs now use quality filtering

### References
- Detailed release notes: `note/v0.14.0-work-BASE_QUALITY_FILTERING.md`
- Performance test: `runs/08-MaP/test_base_quality_performance.sh`
- Test results: `runs/08-MaP/performance_test/`

# [0.13.0] - 2025-01-XX

### Added
- **Distribution-Based K-Factor Prediction**: New advanced k-factor prediction method using statistical distribution analysis
  - New `--k-prediction-method` parameter with options: `background` (default) and `distribution`
  - Uses percentile thresholds and distribution matching to predict optimal k values
  - Implements hierarchical loss function based on feature position in target distribution
  - Automatic fallback to background method if distribution prediction fails
  - Significantly improves k prediction accuracy, especially for 2A3 samples
  - Average k difference reduced from 0.1758 to 0.1175
  - Average AUC gap improved from -1.20% to -1.39%
  - Average vs md_sf improvement: +2.27%

### Changed
- **K-Factor Calculation**: Enhanced reactivity calculation to support multiple k prediction methods
  - Traditional background region method remains default for backward compatibility
  - Distribution method available as opt-in feature
  - Improved logging to indicate which method was used

### Impact
- **Accuracy Improvement**: Distribution method provides better k prediction for many sample types
  - 2A3 human: k difference reduced by 77% (0.354 → 0.080)
  - 2A3 EC: k difference reduced to 0.040 (best performance)
  - Overall: More consistent performance across different sample types
- **Backward Compatible**: Fully backward compatible; default behavior unchanged
  - Existing workflows continue to work without modification
  - New method only activated with explicit `--k-prediction-method distribution` parameter
- **Performance**: Minimal overhead (<1 second) for distribution method
  - Automatic fallback ensures robustness
  - No impact on traditional method performance

### Technical Details
- **File**: `src/reactivity.rs`
- **New Functions**:
  - `predict_k_by_distribution`: Main prediction function using loss minimization
  - `calculate_features_for_k`: Calculates reactivity features for candidate k values
  - `calculate_distribution_loss`: Computes loss based on distribution matching
  - `feature_loss`: Hierarchical loss function for individual features
  - `DistributionTargets`: Target distribution values derived from optimal k analysis
- **Algorithm**: 
  - Uses percentile thresholds (top 10% for high, 40-60% for low mutation rate sites)
  - Searches k range [0.0, 2.0] with 0.01 step
  - Implements weighted loss minimization across three features
- **Target Values**: Based on statistical analysis of optimal k values from 4 diverse samples

### References
- Detailed release notes: `note/v0.13.0-work-DISTRIBUTION_BASED_K_PREDICTION.md`
- Analysis: `runs/08-MaP/distribution_based_method_summary.md`
- Evaluation results: `runs/08-MaP/predicted_k_evaluation_results_distribution_based.csv`

# [0.12.5] - 2025-01-XX

### Fixed
- **Critical Normalization Bug Fix**: Fixed `zarringhalam_remap` function to properly handle negative input values
  - Previous code did not handle negative values, causing negative outputs when using `--linear` parameter
  - Error: Negative values remained negative after remap (e.g., -0.456616 → -0.639262)
  - Now correctly maps all negative values to 0 before applying piecewise linear mapping
  - Ensures all output values are in [0, 1] range as intended by Zarringhalam piecewise linear mapping
  - Fixes normalization for all methods using Siegfried reactivity calculation (`*_sf` methods)

### Changed
- **Zarringhalam Remap Implementation**: Modified `zarringhalam_remap` function to clamp negative values to 0 before mapping
  - Added `x.max(0.0)` operation before piecewise linear mapping
  - All negative input values are now mapped to 0.0
  - Non-negative values maintain original mapping behavior
  - Output values are guaranteed to be in [0, 1] range

### Impact
- **Normalization Output**: All normalized files using `--linear` parameter now have values in [0, 1] range
  - Previously: Output could contain negative values (e.g., -0.639262 to 1.0)
  - Now: All outputs are in [0, 1] range as expected
  - Example: `md_sf_2A3_human_norm.csv` had 926 negative values, now all mapped to [0, 1]
- **Affected Workflows**: 10 Snakefile rules using `--linear` parameter
  - High impact: `rule norm_siegfried` (md_sf) - confirmed affected
  - Medium impact: All `*_sf` Siegfried method rules (4 rules)
  - Low impact: Non-Siegfried method rules (5 rules)
- **Evaluation Metrics**: May change slightly due to corrected normalization
  - Threshold selection may change
  - Sensitivity and specificity balance may shift
  - Requires re-running affected normalization steps
- **Backward Compatible**: Fully backward compatible; no API changes, but data output format improved

### Technical Details
- **File**: `src/norm.rs`
- **Function**: `zarringhalam_remap`
- **Lines Modified**: 170-185
- **Fix Method**: Added `let x_non_neg = x.max(0.0);` before mapping calculation
- **Root Cause**: Function assumed non-negative input, but Siegfried method (v0.12.2) allows negative values

### References
- Bug analysis: `runs/08-MaP/zarringhalam_remap_bug_analysis.md`
- Affected rules: `runs/08-MaP/affected_snakefile_rules.md`
- Detailed release notes: `note/v0.12.5-work-ZARRINGHALAM_REMAP_NEGATIVE_VALUES_FIX.md`

# [0.12.4] - 2025-01-XX

### Fixed
- **Critical Evaluation Bug Fix**: Fixed incorrect sensitivity and specificity calculation in `evaluate_reactivity_accuracy_with_base_matching` function
  - Previous code used fixed threshold `value > 0.0` to calculate confusion matrix, causing incorrect metrics when all reactivity values were positive
  - Error: Sensitivity=1.0, Specificity=0.0 for methods where all reactivity values > 0.0
  - Now uses optimal threshold (threshold with highest F1-score) to calculate all metrics
  - Fixes evaluation metrics for methods: `md`, `br`, `rf`, `rf_norm`, `sm2`, `mpileup` (non-`_sf` versions)

### Changed
- **Evaluation Metrics Calculation**: Modified `evaluate_reactivity_accuracy_with_base_matching` to use optimal threshold selection
  - Calculates F1-score for all possible thresholds
  - Selects threshold with highest F1-score
  - Uses optimal threshold to calculate confusion matrix (TP, FP, TN, FN)
  - All metrics (sensitivity, specificity, accuracy, PPV, NPV, F1) now based on optimal threshold

### Impact
- **Metric Accuracy**: Sensitivity and specificity now correctly reflect classifier performance
  - Example (2A3 EC): md sensitivity improved from 1.0 to 0.7260, specificity from 0.0 to 0.6292
  - Example (2A3 EC): rf_norm sensitivity improved from 1.0 to 0.7713, specificity from 0.0 to 0.6318
- **AUC Values**: Unchanged (AUC considers all thresholds, so it was already correct)
- **F1 Scores**: May change slightly (now more accurate)
- **Backward Compatible**: Fully backward compatible; no API or data format changes

### Technical Details
- **File**: `src/evaluate.rs`
- **Function**: `evaluate_reactivity_accuracy_with_base_matching`
- **Lines Modified**: 784-827
- **Code Logic**: Now consistent with `evaluate_reactivity_accuracy_from_matched_data` function

### References
- Problem analysis: `runs/08-MaP/sensitivity_specificity_issue_analysis.md`
- Fix summary: `runs/08-MaP/evaluation_fix_summary.md`
- Updated results: `runs/08-MaP/auc_summary_simple.md`
- Detailed release notes: `note/v0.12.4-work-EVALUATION_OPTIMAL_THRESHOLD_FIX.md`

# [0.12.3] - 2025-01-XX

### Added
- **Input Format Extension**: Significantly expanded `modtector convert` command's input format support
  - **RNA Framework formats**:
    - `rf-rctools`: Support for RNA Framework / rctools counting output with mutation/RT stop mode distinction via `--rf-count-mutations` parameter
    - `rf-norm`: Support for rf-norm CSV output (reactivity-level data) with validation and NaN preservation
    - `rf-norm-xml`: Support for rf-norm XML output with native Rust XML parsing using `quick_xml` library
  - **ShapeMapper2 format**:
    - `shapemapper-profile`: Support for ShapeMapper2 profile files with dynamic column detection based on header field names (`Nucleotide`, `Sequence`, `Norm_profile`)
    - Robust error handling with skip + count + limited warning strategy for malformed rows
  - **Samtools format**:
    - `samtools-mpileup`: Support for samtools mpileup output (10-column tab-separated format)
    - Automatic mutation count calculation from base counts
    - Strand-specific information support
- **Enhanced Parameter System**: Added/strengthened key parameters to handle format-specific variations
  - `--format/-f`: Manual format specification (with automatic detection as fallback)
  - `--rf-count-mutations`: Distinguish between mutation count and RT stops in rf-count output
  - `--ref-fasta`: Provide reference sequence for formats requiring it (e.g., shapemapper2)
  - `--chr-name`: Specify chromosome/contig name for formats with missing metadata
  - `--strand`: Specify strand orientation when input lacks this information
  - `--log`: Log conversion process warnings and statistics for traceability
- **Format Auto-Detection**: Multi-level detection logic based on file extension, headers, and data patterns
- **Streaming Processing**: Line-by-line processing with progress reporting (every 100k lines) and error tolerance

### Changed
- **Format Detection**: Enhanced `detect_input_format()` function with multi-level detection logic
- **Error Handling**: Improved error handling with error counting, limited warnings (max 10), and continued processing on errors
- **NaN Support**: All conversion functions now properly handle "NaN" string format and preserve NaN semantics

### Technical Details
- **File**: `src/convert.rs`
- **New Functions**: 
  - `parse_rf_rctools_line()`: Parse rf-rctools format
  - `convert_rf_rctools_to_pileup()`: Convert rf-rctools to pileup
  - `convert_rf_norm_to_norm()`: Convert rf-norm CSV to norm
  - `convert_rf_norm_xml_to_norm()`: Parse XML and convert to norm
  - `convert_shapemapper_profile_to_norm()`: Convert shapemapper profile to norm
  - `parse_samtools_mpileup_line()`: Parse samtools mpileup format
  - `convert_samtools_mpileup_to_pileup()`: Convert samtools mpileup to pileup
- **Dependencies**: Added `quick_xml` for XML parsing
- **Code Statistics**: ~500+ lines of new code for format support

### Impact
- **Input Compatibility**: Expanded from single format to 7 supported formats
- **Tool Integration**: Direct integration with RNA Framework, ShapeMapper2, samtools outputs
- **Workflow Simplification**: Reduced external script dependencies, improved maintainability
- **Backward Compatible**: Fully backward compatible with existing functionality

### References
- Detailed release notes: `note/v0.12.3-work-INPUT_FORMAT_EXTENSION.md`

# [0.12.2] - 2025-01-XX

### Changed
- **Critical Accuracy Improvement**: Modified Siegfried method to allow negative reactivity values
  - Removed `.max(0.0)` operation that was forcing all negative values to zero
  - Now matches Shapemapper's implementation: `R_i = T_i - U_i` (allowing negatives)
  - Negative values indicate positions where background mutation rate exceeds modified sample rate
  - This addresses the primary limiting factor for modtector's accuracy

### Impact
- **Data Distribution**: Zero value percentage reduced from 54.8% to ~0.4% (matching Shapemapper)
- **Accuracy**: Expected AUC improvement of ~17.8% (from 0.5795 to ~0.6829)
- **Compatibility**: Maintains backward compatibility; negative values properly handled in downstream steps

### Technical Details
- Updated `src/reactivity.rs` Siegfried method implementation
- Updated documentation in `docs/reactivity_method_description.md`
- Negative values are biologically meaningful and should be preserved

### References
- Analysis documented in `runs/08-MaP/modtector_limitation_analysis_final.md`
- Improvement plan in `runs/08-MaP/modtector_improvement_plan.md`
- Detailed release notes in `note/v0.12.2-work-SIEGFRIED_NEGATIVE_VALUES_SUPPORT.md`

# [0.12.1] - 2025-01-XX

### Fixed
- **Critical Bug Fix**: Fixed NaN value handling in `evaluate_reactivity_accuracy_from_matched_data` function
  - Previous code used `partial_cmp().unwrap_or(Equal)` for sorting, causing panic when NaN values were present
  - Error: "user-provided comparison function does not correctly implement a total order"
  - Now correctly handles NaN values in sorting by placing them at the end
  - Filters out NaN values when calculating total_positive and total_negative counts
  - Skips NaN values when calculating TP/FP/TN/FN metrics
  - Fixes evaluate command failure when processing bamreadcount-derived reactivity files containing NaN values

### Technical Details
- Updated sorting logic in `evaluate.rs` to properly handle NaN values
- Added NaN filtering in statistical calculations
- Ensures consistent behavior with `evaluate_reactivity_accuracy` function which already had proper NaN handling

# [0.12.0] - 2025-01-XX

### Added
- **New `convert` command**: Native Rust implementation for converting bamreadcount format to modtector pileup CSV format
  - Replaces Python script with high-performance Rust module
  - Supports streaming processing for large files (tested up to 35GB)
  - Includes progress reporting and comprehensive error handling
  - Usage: `modtector convert -i input.txt -o output.csv -s +`
- **Bamreadcount reactivity support**: Added `reactivity_bamreadcount` rule in Snakefile workflow
  - Enables reactivity calculation using bamreadcount-derived pileup files
  - Maintains compatibility with existing modtector pileup workflow
  - Both workflows can be used independently

### Changed
- Updated `convert_bamreadcount_to_pileup` Snakefile rule to use native modtector command
  - Removed dependency on external Python script
  - Improved performance and consistency with toolchain

### Technical Details
- New module: `src/convert.rs` with streaming processing implementation
- Uses 8MB buffer for optimal NFS filesystem performance
- Processes ~1-2 million lines/second depending on I/O
- Memory usage: O(1) constant regardless of file size

# [0.11.4] - 2025-12-25

### Fixed
- **Critical Bug Fix**: Fixed SVG template position extraction to prioritize `position.label in template` over `<title>` position
  - SVG templates contain positions in format: `<title>1153 (position.label in template: 1151.A)</title>`
  - Previous code only extracted title position (1153), causing incorrect position mapping
  - Now correctly uses template position (1151) for accurate alignment
- **Critical Bug Fix**: Fixed base-specific shift calculation in SVG position alignment
  - Previous code used global minimum position (often 0) for all bases
  - Now calculates shift using base-specific minimum positions
  - Fixes low matching rates for A and T bases (improved from ~27% to expected >90%)
- Fixed mixed coordinate system handling in SVG templates
  - SVG templates contain both relative positions (0-1870) and genomic coordinates (125931-127798)
  - Genomic coordinates (>= 100000) are now used directly without shift
  - Only relative positions are shifted using base-specific calculations
- Fixed reactivity distribution preview functionality
  - Fixed `color.toHexString is not a function` errors in color picker initialization
  - Improved canvas sizing for high-DPI displays
  - Fixed distribution preview background transparency

### Added
- Reactivity distribution preview graphs above cutoff sliders
  - Histogram-like visualization showing reactivity score distribution
  - Color-coded regions based on cutoff thresholds (low/mid/high)
  - Vertical lines indicating threshold positions
  - Supports both individual base and unified color range modes
  - Transparent background for better integration

### Changed
- Improved position alignment algorithm to handle mixed coordinate systems
- Enhanced error handling in color picker initialization
- Better handling of SVG template position extraction

# [0.11.3] - 2025-12-24

### Added
- SNP filtering in reactivity calculation step with `--snp-cutoff` parameter (default: 0.25)
  - Filters positions with mutation rate >= cutoff before reactivity calculation
  - Filtered positions are marked as NaN in output, preserving position information

### Fixed
- **Critical Bug Fix**: Fixed `chartonum()` function to handle 'U' (uracil) in RNA sequences
  - Previously, RNA reference sequences using U were incorrectly treating all T bases as mutations
  - Changed from `'T' => return 8` to `'T' | 'U' => return 8` to treat U and T as equivalent
  - This fixes 100% mutation rate issues in RNA data where reference uses U but reads contain T
- Fixed sorting panic when NaN values are present in data
  - All sorting functions now correctly handle NaN values (placed at end)
  - Prevents "user-provided comparison function does not correctly implement a total order" panic

### Changed
- NaN value handling throughout the pipeline:
  - Reactivity step: Filtered positions marked as NaN, excluded from k-factor calculation
  - Norm step: NaN values preserved in output, excluded from normalization calculations
  - Evaluate step: NaN values excluded from statistical calculations (AUC, F1-score, etc.)
  - Plot step: NaN values visualized as gray/no data
- All CSV parsing functions now support "NaN" string format
- Statistical calculations (AUC, sensitivity, specificity, etc.) now correctly exclude NaN values

# [0.11.1] - 2025-01-XX

### Fixed
- **Critical Bug Fix**: Corrected stop signal position correction for reverse-strand reads
  - Previously, both forward and reverse strand stop signals were corrected using `start_position - 1`
  - Now correctly applies strand-specific correction:
    - Forward strand: `start_position - 1` (upstream of read start)
    - Reverse strand: `start_position + 1` (downstream of read start)
  - This fix ensures accurate genomic position assignment for RT stop signals according to theoretical principles
  - Affects both `process_region_chunk_ultra_fast` and `process_region_chunk_unified` functions
  - Added boundary checks to ensure corrected positions are within valid reference sequence range

### Changed
- Updated signal correction documentation in `METHOD_DESCRIPTION.md` to accurately reflect strand-specific correction formulas
- Improved code comments to clarify theoretical basis for correction strategy

# [0.11.0] - 2025-01-XX

### Added
- `--batch` parameter for batch processing mode: process multiple BAM files sequentially with glob pattern matching
- `--single-cell` parameter for single-cell unified processing mode: unified processing strategy with automatic cell label extraction
- Glob pattern support for BAM file input (e.g., `/path/to/*sort.bam`) in batch and single-cell modes
- Automatic cell label extraction from BAM filenames (supports RHX pattern and fallback strategies)
- True unified processing: collect all reads from all BAM files for each window, unified pileup with cell label tracking
- Detailed progress reporting: shows chunks processed, processing speed (chunks/s), and estimated time remaining (ETA)
- CPU usage recommendations: automatically suggests optimal CPU count based on BAM file count
- Window-based memory management: process windows sequentially, release memory after each window

### Changed
- Renamed previous `--single-cell` parameter to `--batch` for clarity (it was essentially batch processing)
- Updated `modtector count` to support three modes: standard (single file), batch (sequential), and single-cell unified (parallel)
- Improved output path validation: CSV file for standard mode, directory for batch/single-cell modes
- Enhanced BAM index checking: program now only checks for index existence and prompts user if missing (no automatic indexing)
- Removed data distribution scanning in single-cell mode: directly process all reference sequences without pre-scanning
- Thread count limit: removed hard limit (previously 64), now only shows warning for >128 threads

### Performance
- Single-cell unified mode provides 2-3x performance improvement over batch mode for single-cell data
- Reduced I/O overhead through unified read collection and processing
- Eliminated data distribution scanning overhead (saves 9-19 minutes for 56 files)
- Better memory management: window-based processing prevents memory accumulation
- Improved parallelization efficiency through cross-file unified processing

### Fixed
- Fixed rayon thread pool re-initialization panic in batch processing mode
- Improved error handling for missing BAM index files
- Fixed output path validation for batch/single-cell modes
- Fixed memory accumulation issue in unified processing mode using intermediate files
- Implemented intermediate file mechanism: process chunks → write to temp files → merge and deduplicate → output final CSV

# [0.10.0] - 2025-11-17

### Added
- `modtector duet` sliding-window workflow that fuses normalized stop/mutation reactivity with read-level co-variation (BAM) to infer RNA ensemble compositions without predefining state counts.
- Duet density/window tuning parameters (`--epsilon`, `--min-samples`, `--window-size`, `--window-step`) and enriched window/ensemble summary outputs with co-occurrence statistics.
- Global ensemble aggregation across overlapping windows with new `_global.csv` and `_global_per_base.csv` reports capturing transcript-wide co-variation and per-base read support.
- Confidence scoring for window/global ensembles combining read support, reactivity continuity, and short-range filters to down-weight sparse clusters.
- Multi-threading support for duet workflow with `-t/--threads` parameter (defaults to system CPU count).
- Parallel processing of window analysis, reactivity statistics computation, and global ensemble aggregation.
- Progress reporting for duet workflow: dual-signal reads statistics and high-confidence position assignments.
- Optional `-s/--strand` filter for `modtector count` to limit pileup processing to a single strand when desired.

### Changed
- Renamed the previous `decon` module/command to `duet` (Dynamic Structure Ensemble Decomposition) and refactored outputs to operate at the window level.
- Documentation updated to describe the Duet sliding-window workflow, global ensemble aggregation, new CLI options, and window/global output formats.

### Fixed
- Global ensemble aggregation algorithm: added minimum shared position threshold (5 positions or 15% overlap) to prevent all clusters from being merged into a single ensemble due to transitive connections.
- Improved Union-Find logic in `aggregate_global_ensembles` to require sufficient position overlap before merging clusters.

### Performance
- Duet workflow performance improvements: 3-6x speedup on multi-core systems through parallel window processing.
- Parallelized `analyze_windows`, `compute_window_reactivity_stats`, and `aggregate_global_ensembles` functions using rayon.

## [0.9.7] - 2025-11-12

### Added
- Mod-only reactivity workflow to support smartSHAPE datasets without unmodified controls.
- Command-line and log messaging that clearly indicate when mod-only mode is active.
- Optional `-s/--strand` filter for `modtector count` to limit processing to a single strand when desired.

### Changed
- Documentation examples now describe how to run `modtector reactivity` without `-U`.

### Fixed
- Panic in `modtector norm` when reactivity CSVs contained skipped or malformed rows.


## [0.9.6] - 2025-09-27

### Changed
- Package name changed from `ModDetector` to `modtector`
- Updated all documentation and installation instructions
- Added crates.io installation method

### Added
- System dependency installation guides
- Complete publication preparation checklist

## [0.9.5] - 2025-09-27

### Added
- RNA structure SVG plotting functionality
- Multi-signal support and strand selection
- SVG-only mode for plotting

### Changed
- Enhanced command-line interface with SVG options
- Updated documentation and error handling

### Fixed
- Position mapping and CSV auto-detection improvements

## [0.9.4] - 2025-09-26

### Changed
- Updated help messages and documentation
- Improved code readability and maintainability

## [0.9.3] - 2025-09-25

### Added
- Complete ReadTheDocs documentation
- Comprehensive user guide and examples
- Installation and troubleshooting guides

## [0.9.2] - 2025-09-24

### Changed
- Code cleanup and optimization
- Removed unused functions and dependencies
- Improved error handling

## [0.9.1] - 2025-09-23

### Added
- SVG output improvements
- Enhanced visualization capabilities

## [0.9.0] - 2025-09-22

### Added
- Comprehensive evaluation metrics
- ROC/PR curve generation
- Performance assessment tools

### Fixed
- Evaluation metric calculations
- Accuracy assessment improvements

## [0.8.1] - 2025-09-21

### Added
- Dynamic window scheduling
- Weighted sharding for improved performance

## [0.8.0] - 2025-09-20

### Added
- Optional genome windowing
- Merge functionality improvements

### Changed
- Improved memory management
- Enhanced parallel processing

## [0.7.0] - 2025-09-19

### Changed
- Pileup performance optimizations
- Regression analysis improvements

## [0.6.2] - 2025-09-18

### Added
- DeltaSHAPE implementation
- Advanced comparison methods

## [0.6.1] - 2025-09-17

### Added
- DiffScan analysis features
- Statistical comparison tools

### Fixed
- DiffScan implementation issues

## [0.6.0] - 2025-09-16

### Added
- Compare command implementation
- Sample comparison functionality
- Differential modification detection

## [0.5.9] - 2025-09-15

### Changed
- CLI help completion in English
- Improved user interface

## [0.5.8] - 2025-09-14

### Changed
- Plot simplification
- Unified signal handling

## [0.5.7] - 2025-09-13

### Added
- Enhanced evaluate function
- Improved accuracy assessment

## [0.5.6] - 2025-09-12

### Fixed
- Norm output format issues
- Method combination improvements
- Test script fixes
- Parameter simplification

## [0.5.5] - 2025-09-11

### Changed
- Reactivity output simplification
- Streamlined data processing

## [0.5.4] - 2025-09-10

### Fixed
- Empty strand optimization
- Performance improvements

## [0.5.3] - 2025-09-09

### Added
- Traversal analysis improvements
- Performance monitoring

## [0.5.2] - 2025-09-08

### Added
- Merged table optimization
- Data storage structure analysis

## [0.5.1] - 2025-09-07

### Added
- Reactivity optimization completion
- Performance analysis tools

## [0.5.0] - 2025-09-06

### Added
- Reactivity optimization strategy
- Performance optimization implementation
- Speed optimization analysis

## [0.4.9] - 2025-09-05

### Added
- Logging enhancement
- Improved progress reporting

## [0.4.8] - 2025-09-04

### Added
- English help completion
- Comprehensive documentation

## [0.4.7] - 2025-09-03

### Added
- Parameter validation testing
- Quality assurance improvements

## [0.4.6] - 2025-09-02

### Changed
- Parameter restructure completion
- Improved command-line interface

## [0.4.5] - 2025-09-01

### Added
- Parameter restructure design
- Enhanced usability

## [0.4.4] - 2025-08-31

### Changed
- Parameter naming rationalization
- Improved consistency

## [0.4.3] - 2025-08-30

### Changed
- Removed unused functions
- Code cleanup

## [0.4.2] - 2025-08-29

### Changed
- Simplified reactivity API
- Improved user experience

## [0.4.1] - 2025-08-28

### Changed
- Removed redundant pileup functions
- Streamlined codebase

## [0.4.0] - 2025-08-27

### Added
- Redundancy analysis and planning
- Code optimization strategy

## [0.3.4] - 2025-08-26

### Added
- Progress display fixes
- Evaluate output optimization
- Plot reactivity support
- Output optimization across all modules

### Fixed
- Reactivity method fixes
- AUC output optimization
- Time format unification

## [0.3.3] - 2025-08-25

### Added
- Float issue audit
- Reactivity integration testing
- IO optimization

### Fixed
- Norm float issue resolution

## [0.3.2] - 2025-08-24

### Added
- IO optimization final report
- Parallel pileup implementation

### Fixed
- Optimization summary
- Fixed parallel pileup issues

## [0.3.1] - 2025-08-23

### Added
- Final summary
- Eukaryotic transcript extraction

## [0.3.0] - 2025-08-22

### Added
- GTF support implementation
- GTF/GFF3 comparison analysis

## [0.2.0] - 2025-08-21

### Added
- AUC analysis
- Method combination testing
- Reactivity methods implementation

### Fixed
- Ding method fixes and comparison

## [0.1.0] - 2025-08-20

### Added
- Initial release
- Basic modification detection
- Pileup analysis
- Reactivity calculation
- Basic plotting functionality

### Fixed
- Position correction details
- BAM reads direction analysis
