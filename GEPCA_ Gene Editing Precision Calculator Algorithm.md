# GEPCA: Gene Editing Precision Calculator Algorithm

## Overview

GEPCA (Gene Editing Precision Calculator Algorithm) is a computational framework designed to predict and quantify the precision of CRISPR-Cas9 gene editing outcomes. By integrating sequence analysis, structural prediction, machine learning, and statistical modeling, GEPCA provides researchers with a comprehensive tool for optimizing gene editing experiments, minimizing off-target effects, and maximizing editing efficiency.

## Installation

```bash
# Clone the repository
git clone https://github.com/gepca/gepca.git

# Navigate to the directory
cd gepca

# Install the package
pip install -e .
```

## Quick Start

```python
# Import the GEPCA modules
from gepca import GeneEditingPrecisionCalculator
from gepca.visualizer import GEPCAVisualizer

# Initialize the calculator with the reference genome
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")

# Define a target sequence and guide RNA
target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
guide_rna = "ACGTACGTACGTACGTACGT"

# Perform analysis
results = calculator.analyze(target_sequence, guide_rna)

# Print results
print(f"Precision Score: {results['precision_score']['precision_score']:.4f}")
print(f"Confidence Interval: {results['precision_score']['confidence_interval']}")
print(f"On-Target Efficiency: {results['editing_outcomes']['on_target_efficiency']:.4f}")
print(f"Off-Target Probability: {results['off_target_analysis']['total_probability']:.4f}")

# Visualize results
visualizer = GEPCAVisualizer(theme='default')
visualizer.plot_precision_score(results)
```

## Core Components

GEPCA consists of four main components:

### 1. Target Sequence Analyzer

The Target Sequence Analyzer evaluates the primary target sequence and identifies potential off-target sites based on sequence homology. Key features include:

- **PAM Site Identification**: Locates all PAM sites (NGG for SpCas9) within the target genome
- **Homology Search**: Identifies sequences with varying degrees of homology to the target
- **Seed Region Analysis**: Gives special consideration to the seed region (8-12 bp proximal to PAM)
- **Repetitive Element Detection**: Flags repetitive elements that may lead to multiple off-target effects

### 2. Structural Context Evaluator

The Structural Context Evaluator assesses the accessibility of the target site and potential off-target sites. Key features include:

- **Chromatin State Assessment**: Evaluates chromatin accessibility at target and off-target sites
- **DNA Secondary Structure Prediction**: Predicts DNA secondary structures that may affect binding
- **Epigenetic Modification Analysis**: Considers the impact of DNA methylation and histone modifications
- **Local Sequence Context**: Analyzes the impact of surrounding sequences on editing efficiency

### 3. Deep Learning Prediction Engine

The Deep Learning Prediction Engine utilizes state-of-the-art neural network architectures to predict editing outcomes. Key features include:

- **On-Target Efficiency Prediction**: Predicts the efficiency of editing at the target site
- **Off-Target Probability Calculation**: Calculates the probability of editing at potential off-target sites
- **Edit Pattern Prediction**: Predicts the types of edits (insertions, deletions, substitutions) likely to occur
- **Variant-Specific Modeling**: Incorporates models specific to different Cas9 variants and other nucleases

### 4. Statistical Confidence Calculator

The Statistical Confidence Calculator provides confidence scores and uncertainty estimates for predictions. Key features include:

- **Confidence Interval Generation**: Generates confidence intervals for efficiency predictions
- **Uncertainty Quantification**: Quantifies uncertainty in off-target predictions
- **Sensitivity Analysis**: Assesses the sensitivity of predictions to input parameters
- **Comparative Benchmarking**: Compares predictions to experimental data for validation

## API Reference

### GeneEditingPrecisionCalculator

The main class for performing GEPCA analysis.

#### Constructor

```python
GeneEditingPrecisionCalculator(reference_genome, weights=None)
```

- `reference_genome`: Path to the reference genome file
- `weights`: Optional dictionary of weights for the precision score calculation. Default weights are:
  - `on_target_efficiency`: 0.4
  - `off_target_probability`: 0.3
  - `structural_context`: 0.15
  - `edit_characteristics`: 0.15

#### Methods

##### analyze

```python
analyze(target_sequence, guide_rna, cell_type='HEK293', cas_variant='SpCas9')
```

Perform a complete analysis of the gene editing precision.

- `target_sequence`: The target DNA sequence
- `guide_rna`: The guide RNA sequence
- `cell_type`: The cell type for context-specific evaluation (default: 'HEK293')
- `cas_variant`: The Cas variant being used (default: 'SpCas9')

Returns a dictionary containing comprehensive analysis results.

##### analyze_target_sequence

```python
analyze_target_sequence(target_sequence, guide_rna)
```

Analyze the target sequence and identify potential off-target sites.

- `target_sequence`: The DNA sequence containing the target site
- `guide_rna`: The guide RNA sequence (without PAM)

Returns a dictionary containing analysis results.

##### evaluate_structural_context

```python
evaluate_structural_context(target_site, cell_type)
```

Evaluate the structural context of the target site.

- `target_site`: Dictionary containing target site information
- `cell_type`: The cell type for context-specific evaluation

Returns a dictionary containing structural context evaluation.

##### predict_editing_outcomes

```python
predict_editing_outcomes(guide_rna, target_sequence, cas_variant='SpCas9')
```

Predict the editing outcomes using deep learning models.

- `guide_rna`: The guide RNA sequence
- `target_sequence`: The target DNA sequence
- `cas_variant`: The Cas variant being used (default: 'SpCas9')

Returns a dictionary containing editing outcome predictions.

##### calculate_off_target_probability

```python
calculate_off_target_probability(guide_rna, off_target_sites)
```

Calculate the probability of off-target effects.

- `guide_rna`: The guide RNA sequence
- `off_target_sites`: List of potential off-target sites

Returns a dictionary containing off-target probabilities.

##### calculate_precision_score

```python
calculate_precision_score(on_target_efficiency, off_target_probability, structural_score, edit_char_score)
```

Calculate the overall precision score.

- `on_target_efficiency`: Predicted on-target efficiency
- `off_target_probability`: Predicted off-target probability
- `structural_score`: Structural context score
- `edit_char_score`: Edit characteristics score

Returns a dictionary containing precision score and components.

### GEPCAVisualizer

The visualization class for GEPCA results.

#### Constructor

```python
GEPCAVisualizer(theme='default')
```

- `theme`: Visual theme for plots ('default', 'dark', or 'publication')

#### Methods

##### plot_precision_score

```python
plot_precision_score(results, save_path=None)
```

Plot the precision score and its components.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

##### plot_off_target_analysis

```python
plot_off_target_analysis(results, save_path=None)
```

Plot the off-target analysis results.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

##### plot_edit_patterns

```python
plot_edit_patterns(results, save_path=None)
```

Plot the predicted edit patterns.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

##### plot_target_sequence_analysis

```python
plot_target_sequence_analysis(results, save_path=None)
```

Plot the target sequence analysis.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

##### plot_structural_context

```python
plot_structural_context(results, save_path=None)
```

Plot the structural context analysis.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

##### create_summary_dashboard

```python
create_summary_dashboard(results, save_path=None)
```

Create a comprehensive dashboard of all analysis results.

- `results`: Results dictionary from the GEPCA analysis
- `save_path`: Optional path to save the figure

Returns a matplotlib figure.

## Examples

### Basic Usage

```python
from gepca import GeneEditingPrecisionCalculator
from gepca.visualizer import GEPCAVisualizer

# Initialize the calculator
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")

# Define a target sequence and guide RNA
target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
guide_rna = "ACGTACGTACGTACGTACGT"

# Perform analysis
results = calculator.analyze(target_sequence, guide_rna)

# Print results
print(f"Precision Score: {results['precision_score']['precision_score']:.4f}")
print(f"Confidence Interval: {results['precision_score']['confidence_interval']}")
print(f"On-Target Efficiency: {results['editing_outcomes']['on_target_efficiency']:.4f}")
print(f"Off-Target Probability: {results['off_target_analysis']['total_probability']:.4f}")

# Initialize the visualizer
visualizer = GEPCAVisualizer(theme='default')

# Create and save visualizations
visualizer.plot_precision_score(results, save_path="precision_score.png")
visualizer.plot_off_target_analysis(results, save_path="off_target_analysis.png")
visualizer.plot_edit_patterns(results, save_path="edit_patterns.png")
visualizer.plot_target_sequence_analysis(results, save_path="target_sequence.png")
visualizer.plot_structural_context(results, save_path="structural_context.png")

# Create comprehensive dashboard
visualizer.create_summary_dashboard(results, save_path="gepca_dashboard.png")
```

### Comparing Multiple Guide RNAs

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gepca import GeneEditingPrecisionCalculator

# Initialize the calculator
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")

# Define multiple guide RNAs targeting the same region
guide_rnas = [
    "ACGTACGTACGTACGTACGT",
    "ACGTACGTACGTACGTTCGT",
    "ACGTACGTACGTTCGTACGT",
    "TCGTACGTACGTACGTACGT",
    "ACGTTCGTACGTACGTACGT"
]

# Define target sequence
target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

# Compare guide RNAs
comparison_results = []
for guide_rna in guide_rnas:
    results = calculator.analyze(target_sequence, guide_rna)
    comparison_results.append({
        'Guide RNA': guide_rna,
        'Precision Score': round(results['precision_score']['precision_score'], 3),
        'On-Target Efficiency': round(results['editing_outcomes']['on_target_efficiency'], 3),
        'Off-Target Probability': round(results['off_target_analysis']['total_probability'], 3),
        'Confidence Interval': (round(results['precision_score']['confidence_interval'][0], 3),
                               round(results['precision_score']['confidence_interval'][1], 3))
    })

# Create a DataFrame for easy comparison
comparison_df = pd.DataFrame(comparison_results)
print(comparison_df)

# Visualize comparison
plt.figure(figsize=(10, 6))

# Create bar chart for precision scores
x = np.arange(len(guide_rnas))
width = 0.25

plt.bar(x - width, comparison_df['Precision Score'], width, label='Precision Score', color='#4285F4')
plt.bar(x, comparison_df['On-Target Efficiency'], width, label='On-Target Efficiency', color='#34A853')
plt.bar(x + width, comparison_df['Off-Target Probability'], width, label='Off-Target Probability', color='#EA4335')

plt.xlabel('Guide RNA')
plt.ylabel('Score')
plt.title('Comparison of Guide RNAs')
plt.xticks(x, [f'Guide {i+1}' for i in range(len(guide_rnas))])
plt.legend()
plt.ylim(0, 1)

plt.tight_layout()
plt.savefig('guide_rna_comparison.png', dpi=300)
plt.show()
```

### Cell-Type Specific Analysis

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gepca import GeneEditingPrecisionCalculator

# Initialize the calculator
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")

# Define cell types
cell_types = ['HEK293', 'K562', 'HepG2', 'iPSC', 'CD34+']

# Define target sequence and guide RNA
target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
guide_rna = "ACGTACGTACGTACGTACGT"

# Perform cell-type specific analysis
cell_type_results = []
for cell_type in cell_types:
    results = calculator.analyze(target_sequence, guide_rna, cell_type=cell_type)
    cell_type_results.append({
        'Cell Type': cell_type,
        'Precision Score': round(results['precision_score']['precision_score'], 3),
        'On-Target Efficiency': round(results['editing_outcomes']['on_target_efficiency'], 3),
        'Structural Context Score': round(results['structural_context']['structural_score'], 3)
    })

# Create a DataFrame for easy comparison
cell_type_df = pd.DataFrame(cell_type_results)
print(cell_type_df)

# Visualize cell-type specific analysis
plt.figure(figsize=(10, 6))

# Create bar chart for cell-type specific results
x = np.arange(len(cell_types))
width = 0.25

plt.bar(x - width, cell_type_df['Precision Score'], width, label='Precision Score', color='#4285F4')
plt.bar(x, cell_type_df['On-Target Efficiency'], width, label='On-Target Efficiency', color='#34A853')
plt.bar(x + width, cell_type_df['Structural Context Score'], width, label='Structural Context Score', color='#FBBC05')

plt.xlabel('Cell Type')
plt.ylabel('Score')
plt.title('Cell-Type Specific Analysis')
plt.xticks(x, cell_types)
plt.legend()
plt.ylim(0, 1)

plt.tight_layout()
plt.savefig('cell_type_analysis.png', dpi=300)
plt.show()
```

## Advanced Usage

### Custom Weighting Factors

You can customize the weighting factors used in the precision score calculation:

```python
# Define custom weights
custom_weights = {
    'on_target_efficiency': 0.5,    # Increase weight for on-target efficiency
    'off_target_probability': 0.3,  # Keep the same weight for off-target probability
    'structural_context': 0.1,      # Decrease weight for structural context
    'edit_characteristics': 0.1     # Decrease weight for edit characteristics
}

# Initialize the calculator with custom weights
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa", weights=custom_weights)

# Perform analysis with custom weights
results = calculator.analyze(target_sequence, guide_rna)
```

### Batch Processing

For analyzing multiple guide RNAs in batch:

```python
import pandas as pd
from gepca import GeneEditingPrecisionCalculator

# Initialize the calculator
calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")

# Load guide RNAs from a file
guide_rnas_df = pd.read_csv('guide_rnas.csv')

# Process each guide RNA
results_list = []
for index, row in guide_rnas_df.iterrows():
    guide_rna = row['guide_rna']
    target_sequence = row['target_sequence']
    
    # Perform analysis
    results = calculator.analyze(target_sequence, guide_rna)
    
    # Extract key metrics
    results_list.append({
        'guide_rna': guide_rna,
        'precision_score': results['precision_score']['precision_score'],
        'on_target_efficiency': results['editing_outcomes']['on_target_efficiency'],
        'off_target_probability': results['off_target_analysis']['total_probability']
    })

# Create a DataFrame with results
results_df = pd.DataFrame(results_list)

# Save results to a CSV file
results_df.to_csv('gepca_results.csv', index=False)
```

## Troubleshooting

### Common Issues

1. **Missing Reference Genome**
   
   Error: `FileNotFoundError: [Errno 2] No such file or directory: 'hg38.fa'`
   
   Solution: Ensure the reference genome file is in the correct location or provide the full path to the file.

2. **Invalid Guide RNA**
   
   Error: `ValueError: Invalid guide RNA sequence`
   
   Solution: Ensure the guide RNA is a valid 20-nucleotide sequence containing only A, T, G, and C.

3. **Memory Issues**
   
   Error: `MemoryError`
   
   Solution: For large-scale analyses, consider processing guide RNAs in smaller batches or using a machine with more memory.

### Getting Help

If you encounter any issues or have questions, please:

1. Check the [documentation](https://gepca.readthedocs.io/)
2. Look for similar issues in the [GitHub repository](https://github.com/gepca/gepca/issues)
3. Open a new issue if your problem is not already addressed

## Citation

If you use GEPCA in your research, please cite:

```
Author, A. (2025). Gene Editing Precision Calculator Algorithm: A Novel Computational Framework for Predicting CRISPR-Cas9 Editing Outcomes. Journal of Computational Biology, XX(X), XXX-XXX.
```

## License

GEPCA is released under the MIT License. See the LICENSE file for details.

