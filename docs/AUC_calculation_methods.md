# AUC Calculation Methods

## Problem Description

In ModDetector's evaluation functionality, there may be differences between the AUC values in statistical results and the AUC values displayed in plots. This is because the code implements two different AUC calculation methods.

## Two AUC Calculation Methods

### Method 1: Simple Sorting Method

**Function Location**: `src/evaluate.rs` lines 707-713  
**Used by Functions**: `evaluate_reactivity_accuracy_with_auto_shift` and `evaluate_reactivity_accuracy_with_base_matching`

**Calculation Steps**:
1. Calculate ROC points (FPR, TPR) for each unique reactivity value
2. Sort ROC points by FPR
3. Directly calculate AUC using trapezoidal rule:
   ```rust
   auc += (fpr_curr - fpr_prev) * (tpr_curr + tpr_prev) / 2.0;
   ```

**Characteristics**:
- Simple and direct
- May not correctly handle cases where there are identical FPR values but different TPR values
- May result in non-monotonic ROC curves (TPR may decrease)

### Method 2: Monotonic ROC Points Method

**Function Location**: `src/evaluate.rs` lines 1851-1890  
**Used by Function**: `evaluate_reactivity_accuracy_from_matched_data`

**Calculation Steps**:
1. Calculate ROC points (FPR, TPR) for each unique reactivity value
2. Sort ROC points by FPR
3. **Key Step**: Ensure monotonicity of ROC points
   - For identical FPR values, take the maximum TPR value
   - Create a monotonic ROC point sequence
4. Calculate AUC using monotonic ROC points:
   ```rust
   auc += (fpr_curr - fpr_prev) * (tpr_prev + tpr_curr) / 2.0;
   ```

**Characteristics**:
- Ensures ROC curve is monotonically increasing (TPR does not decrease)
- Conforms to standard ROC curve definition
- More accurately reflects classifier performance

## Actual Usage

### Overall Evaluation

When using `evaluate_both_signals_main` or `evaluate_reactivity_accuracy_with_auto_shift_main`:
- Uses **Method 1** (Simple Sorting Method)
- Generated files: `evaluation_stop.txt`, `evaluation_mutation.txt`
- Displayed AUC values: e.g., `0.5141` (stop), `0.5415` (mutation)

### Base-Specific Evaluation

When using `evaluate_reactivity_accuracy_by_base`:
- Uses **Method 2** (Monotonic ROC Points Method)
- Generated files: `evaluation_by_base_*_*.overall.txt`
- Displayed AUC values: e.g., `0.484` (stop overall)

## Why Are There Differences?

1. **Different Handling of Identical FPR Values**:
   - Method 1: May retain multiple points with identical FPR but different TPR
   - Method 2: Takes maximum TPR for identical FPR, ensuring monotonicity

2. **ROC Curve Definition**:
   - Standard ROC curves should be monotonic (TPR increases or remains constant as FPR increases)
   - Method 2 better conforms to this standard

3. **Numerical Precision**:
   - Due to different handling approaches, the final calculated AUC values will have slight differences

## Recommendations

1. **Unify to Use Monotonic ROC Points Method**:
   - This is a more standard and accurate method
   - Recommend modifying the `evaluate_reactivity_accuracy_with_auto_shift` function to also use the monotonic ROC points method

2. **Documentation**:
   - Clearly annotate which AUC calculation method is used by different evaluation methods
   - Help users understand result differences

3. **Code Refactoring**:
   - Extract AUC calculation logic as an independent function
   - Unify to use monotonic ROC points method

## Example

Assume the following ROC points:
- (0.1, 0.2)
- (0.1, 0.3)  ← Same FPR, different TPR
- (0.2, 0.4)

**Method 1**: Retains all points, may calculate as:
```
AUC = (0.1-0) * (0.2+0)/2 + (0.1-0.1) * (0.3+0.2)/2 + (0.2-0.1) * (0.4+0.3)/2
```

**Method 2**: Takes maximum TPR for identical FPR, calculates as:
```
Monotonic points: (0.1, 0.3), (0.2, 0.4)
AUC = (0.1-0) * (0.3+0)/2 + (0.2-0.1) * (0.4+0.3)/2
```

## Related Code Locations

- `evaluate_reactivity_accuracy_with_auto_shift`: lines 1483-1489
- `evaluate_reactivity_accuracy_with_base_matching`: lines 707-713
- `evaluate_reactivity_accuracy_from_matched_data`: lines 1851-1890
