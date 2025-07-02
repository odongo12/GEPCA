"""
Gene Editing Precision Calculator Algorithm (GEPCA)
A computational framework for predicting and quantifying CRISPR-Cas9 gene editing precision.

This implementation demonstrates the core components and workflow of the algorithm.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import re

class GeneEditingPrecisionCalculator:
    """
    Main class for the Gene Editing Precision Calculator Algorithm.
    """
    
    def __init__(self, reference_genome: str, weights: Optional[Dict[str, float]] = None):
        """
        Initialize the Gene Editing Precision Calculator.
        
        Args:
            reference_genome: Path to the reference genome file
            weights: Optional dictionary of weights for the precision score calculation
        """
        self.reference_genome = reference_genome
        self.weights = weights or {
            'on_target_efficiency': 0.4,
            'off_target_probability': 0.3,
            'structural_context': 0.15,
            'edit_characteristics': 0.15
        }
        
        # Load pre-trained models (in a real implementation)
        self.on_target_model = self._load_on_target_model()
        self.off_target_model = self._load_off_target_model()
        self.structural_model = self._load_structural_model()
        self.edit_pattern_model = self._load_edit_pattern_model()
        
    def _load_on_target_model(self):
        """Load the pre-trained on-target efficiency prediction model."""
        # In a real implementation, this would load a trained deep learning model
        # For demonstration, we'll use a placeholder function
        return lambda x: np.random.uniform(0.5, 0.9)  # Placeholder
        
    def _load_off_target_model(self):
        """Load the pre-trained off-target prediction model."""
        # In a real implementation, this would load a trained deep learning model
        return lambda x, y: np.random.uniform(0.0, 0.5)  # Placeholder
        
    def _load_structural_model(self):
        """Load the pre-trained structural context prediction model."""
        # In a real implementation, this would load a trained model
        return lambda x: np.random.uniform(0.3, 0.8)  # Placeholder
        
    def _load_edit_pattern_model(self):
        """Load the pre-trained edit pattern prediction model."""
        # In a real implementation, this would load a trained model
        return lambda x: {
            'insertion': np.random.uniform(0.2, 0.4),
            'deletion': np.random.uniform(0.4, 0.6),
            'substitution': np.random.uniform(0.1, 0.2)
        }  # Placeholder
    
    def analyze_target_sequence(self, target_sequence: str, guide_rna: str) -> Dict:
        """
        Analyze the target sequence and identify potential off-target sites.
        
        Args:
            target_sequence: The DNA sequence containing the target site
            guide_rna: The guide RNA sequence (without PAM)
            
        Returns:
            Dictionary containing analysis results
        """
        # Check if guide RNA is valid
        if not self._validate_guide_rna(guide_rna):
            raise ValueError("Invalid guide RNA sequence")
        
        # Find PAM sites in the target sequence
        pam_sites = self._find_pam_sites(target_sequence)
        
        # Identify the target site
        target_site = self._identify_target_site(target_sequence, guide_rna)
        
        # Find potential off-target sites
        off_target_sites = self._find_potential_off_targets(guide_rna)
        
        # Analyze seed region
        seed_region_analysis = self._analyze_seed_region(guide_rna)
        
        return {
            'pam_sites': pam_sites,
            'target_site': target_site,
            'off_target_sites': off_target_sites,
            'seed_region_analysis': seed_region_analysis
        }
    
    def _validate_guide_rna(self, guide_rna: str) -> bool:
        """
        Validate the guide RNA sequence.
        
        Args:
            guide_rna: The guide RNA sequence
            
        Returns:
            Boolean indicating if the guide RNA is valid
        """
        # Check length (typically 20 nucleotides for SpCas9)
        if len(guide_rna) != 20:
            return False
        
        # Check if it contains only valid nucleotides
        if not all(n in 'ATGC' for n in guide_rna.upper()):
            return False
        
        return True
    
    def _find_pam_sites(self, sequence: str) -> List[int]:
        """
        Find all PAM sites (NGG for SpCas9) in the sequence.
        
        Args:
            sequence: The DNA sequence to search
            
        Returns:
            List of positions where PAM sites are found
        """
        pam_pattern = re.compile(r'[ATGC]GG')
        return [match.start() for match in pam_pattern.finditer(sequence.upper())]
    
    def _identify_target_site(self, sequence: str, guide_rna: str) -> Dict:
        """
        Identify the target site in the sequence.
        
        Args:
            sequence: The DNA sequence containing the target site
            guide_rna: The guide RNA sequence
            
        Returns:
            Dictionary with target site information
        """
        # In a real implementation, this would perform a more sophisticated search
        # For demonstration, we'll assume the guide RNA perfectly matches somewhere in the sequence
        guide_rna = guide_rna.upper()
        sequence = sequence.upper()
        
        position = sequence.find(guide_rna)
        if position == -1:
            # Try the reverse complement
            rev_comp = self._reverse_complement(guide_rna)
            position = sequence.find(rev_comp)
            if position == -1:
                return {'found': False}
            else:
                return {
                    'found': True,
                    'position': position,
                    'strand': 'reverse',
                    'sequence': rev_comp
                }
        else:
            return {
                'found': True,
                'position': position,
                'strand': 'forward',
                'sequence': guide_rna
            }
    
    def _reverse_complement(self, sequence: str) -> str:
        """
        Get the reverse complement of a DNA sequence.
        
        Args:
            sequence: The DNA sequence
            
        Returns:
            The reverse complement sequence
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(sequence.upper()))
    
    def _find_potential_off_targets(self, guide_rna: str) -> List[Dict]:
        """
        Find potential off-target sites for the guide RNA.
        
        Args:
            guide_rna: The guide RNA sequence
            
        Returns:
            List of dictionaries containing off-target site information
        """
        # In a real implementation, this would search the reference genome
        # For demonstration, we'll generate some synthetic off-target sites
        off_targets = []
        
        # Generate 5 synthetic off-target sites with varying mismatches
        for i in range(5):
            mismatches = np.random.randint(1, 4)
            mismatch_positions = np.random.choice(range(len(guide_rna)), mismatches, replace=False)
            
            off_target_seq = list(guide_rna.upper())
            for pos in mismatch_positions:
                current_base = off_target_seq[pos]
                # Replace with a different base
                options = [b for b in 'ATGC' if b != current_base]
                off_target_seq[pos] = np.random.choice(options)
            
            off_targets.append({
                'sequence': ''.join(off_target_seq),
                'mismatches': mismatches,
                'mismatch_positions': sorted(mismatch_positions),
                'chromosome': f'chr{np.random.randint(1, 23)}',
                'position': np.random.randint(1, 1000000)
            })
        
        return off_targets
    
    def _analyze_seed_region(self, guide_rna: str) -> Dict:
        """
        Analyze the seed region of the guide RNA.
        
        Args:
            guide_rna: The guide RNA sequence
            
        Returns:
            Dictionary containing seed region analysis
        """
        # The seed region is typically the 8-12 bp proximal to the PAM site
        seed_region = guide_rna[-12:-1].upper()
        
        # Calculate GC content
        gc_content = (seed_region.count('G') + seed_region.count('C')) / len(seed_region)
        
        return {
            'seed_sequence': seed_region,
            'gc_content': gc_content,
            'length': len(seed_region)
        }
    
    def evaluate_structural_context(self, target_site: Dict, cell_type: str) -> Dict:
        """
        Evaluate the structural context of the target site.
        
        Args:
            target_site: Dictionary containing target site information
            cell_type: The cell type for context-specific evaluation
            
        Returns:
            Dictionary containing structural context evaluation
        """
        # In a real implementation, this would use epigenomic data for the specific cell type
        # For demonstration, we'll generate synthetic scores
        
        chromatin_accessibility = np.random.uniform(0.5, 1.0)
        dna_shape_score = np.random.uniform(0.3, 0.9)
        methylation_status = np.random.choice(['methylated', 'unmethylated'], 
                                             p=[0.3, 0.7])
        
        # Calculate overall structural context score
        structural_score = 0.5 * chromatin_accessibility + 0.5 * dna_shape_score
        if methylation_status == 'methylated':
            structural_score *= 0.8  # Penalty for methylation
        
        return {
            'chromatin_accessibility': chromatin_accessibility,
            'dna_shape_score': dna_shape_score,
            'methylation_status': methylation_status,
            'structural_score': structural_score
        }
    
    def predict_editing_outcomes(self, guide_rna: str, target_sequence: str, 
                                cas_variant: str = 'SpCas9') -> Dict:
        """
        Predict the editing outcomes using deep learning models.
        
        Args:
            guide_rna: The guide RNA sequence
            target_sequence: The target DNA sequence
            cas_variant: The Cas variant being used
            
        Returns:
            Dictionary containing editing outcome predictions
        """
        # Predict on-target efficiency
        on_target_efficiency = self.on_target_model(guide_rna)
        
        # Predict edit patterns
        edit_patterns = self.edit_pattern_model(guide_rna)
        
        # Calculate edit characteristics score
        desired_pattern = 'deletion'  # This would be specified by the user in a real implementation
        edit_char_score = self._calculate_edit_characteristics_score(edit_patterns, desired_pattern)
        
        return {
            'on_target_efficiency': on_target_efficiency,
            'edit_patterns': edit_patterns,
            'edit_characteristics_score': edit_char_score
        }
    
    def _calculate_edit_characteristics_score(self, edit_patterns: Dict, desired_pattern: str) -> float:
        """
        Calculate the edit characteristics score based on the predicted edit patterns.
        
        Args:
            edit_patterns: Dictionary of predicted edit pattern probabilities
            desired_pattern: The desired edit pattern
            
        Returns:
            Edit characteristics score
        """
        # Simple implementation: score is the probability of the desired pattern
        return edit_patterns.get(desired_pattern, 0.0)
    
    def calculate_off_target_probability(self, guide_rna: str, off_target_sites: List[Dict]) -> Dict:
        """
        Calculate the probability of off-target effects.
        
        Args:
            guide_rna: The guide RNA sequence
            off_target_sites: List of potential off-target sites
            
        Returns:
            Dictionary containing off-target probabilities
        """
        site_probabilities = {}
        total_probability = 0.0
        
        for site in off_target_sites:
            # Calculate probability for each site
            prob = self.off_target_model(guide_rna, site['sequence'])
            
            # Apply mismatch penalty based on number and position of mismatches
            mismatch_penalty = self._calculate_mismatch_penalty(site['mismatch_positions'])
            prob *= mismatch_penalty
            
            site_id = f"{site['chromosome']}:{site['position']}"
            site_probabilities[site_id] = prob
            total_probability += prob
        
        # Normalize if total probability exceeds 1
        if total_probability > 1.0:
            for site_id in site_probabilities:
                site_probabilities[site_id] /= total_probability
            total_probability = 1.0
        
        return {
            'site_probabilities': site_probabilities,
            'total_probability': total_probability
        }
    
    def _calculate_mismatch_penalty(self, mismatch_positions: List[int]) -> float:
        """
        Calculate a penalty factor based on mismatch positions.
        
        Args:
            mismatch_positions: List of positions where mismatches occur
            
        Returns:
            Penalty factor (0-1)
        """
        # Mismatches in the seed region (positions 12-20) are more tolerated
        # than mismatches in the non-seed region (positions 1-11)
        seed_mismatches = sum(1 for pos in mismatch_positions if pos >= 12)
        non_seed_mismatches = len(mismatch_positions) - seed_mismatches
        
        # Calculate penalty (higher penalty = lower probability)
        seed_penalty = 0.7 ** seed_mismatches
        non_seed_penalty = 0.3 ** non_seed_mismatches
        
        return seed_penalty * non_seed_penalty
    
    def calculate_precision_score(self, on_target_efficiency: float, off_target_probability: float,
                                 structural_score: float, edit_char_score: float) -> Dict:
        """
        Calculate the overall precision score.
        
        Args:
            on_target_efficiency: Predicted on-target efficiency
            off_target_probability: Predicted off-target probability
            structural_score: Structural context score
            edit_char_score: Edit characteristics score
            
        Returns:
            Dictionary containing precision score and components
        """
        # Calculate weighted precision score
        precision_score = (
            on_target_efficiency * self.weights['on_target_efficiency'] -
            off_target_probability * self.weights['off_target_probability'] +
            structural_score * self.weights['structural_context'] +
            edit_char_score * self.weights['edit_characteristics']
        )
        
        # Ensure score is between 0 and 1
        precision_score = max(0.0, min(1.0, precision_score))
        
        # Calculate confidence interval (simplified)
        confidence_interval = (
            max(0.0, precision_score - 0.1),
            min(1.0, precision_score + 0.1)
        )
        
        return {
            'precision_score': precision_score,
            'confidence_interval': confidence_interval,
            'components': {
                'on_target_efficiency': on_target_efficiency,
                'off_target_probability': off_target_probability,
                'structural_score': structural_score,
                'edit_char_score': edit_char_score
            }
        }
    
    def analyze(self, target_sequence: str, guide_rna: str, cell_type: str = 'HEK293',
               cas_variant: str = 'SpCas9') -> Dict:
        """
        Perform a complete analysis of the gene editing precision.
        
        Args:
            target_sequence: The target DNA sequence
            guide_rna: The guide RNA sequence
            cell_type: The cell type for context-specific evaluation
            cas_variant: The Cas variant being used
            
        Returns:
            Dictionary containing comprehensive analysis results
        """
        # Step 1: Analyze target sequence
        sequence_analysis = self.analyze_target_sequence(target_sequence, guide_rna)
        
        # Step 2: Evaluate structural context
        structural_context = self.evaluate_structural_context(
            sequence_analysis['target_site'], cell_type)
        
        # Step 3: Predict editing outcomes
        editing_outcomes = self.predict_editing_outcomes(guide_rna, target_sequence, cas_variant)
        
        # Step 4: Calculate off-target probability
        off_target_analysis = self.calculate_off_target_probability(
            guide_rna, sequence_analysis['off_target_sites'])
        
        # Step 5: Calculate precision score
        precision_score = self.calculate_precision_score(
            editing_outcomes['on_target_efficiency'],
            off_target_analysis['total_probability'],
            structural_context['structural_score'],
            editing_outcomes['edit_characteristics_score']
        )
        
        # Compile results
        return {
            'sequence_analysis': sequence_analysis,
            'structural_context': structural_context,
            'editing_outcomes': editing_outcomes,
            'off_target_analysis': off_target_analysis,
            'precision_score': precision_score
        }


# Example usage
if __name__ == "__main__":
    # Initialize the calculator
    calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")
    
    # Example target sequence and guide RNA
    target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    guide_rna = "ACGTACGTACGTACGTACGT"
    
    # Perform analysis
    results = calculator.analyze(target_sequence, guide_rna)
    
    # Print results
    print(f"Precision Score: {results['precision_score']['precision_score']:.4f}")
    print(f"Confidence Interval: {results['precision_score']['confidence_interval']}")
    print(f"On-Target Efficiency: {results['editing_outcomes']['on_target_efficiency']:.4f}")
    print(f"Off-Target Probability: {results['off_target_analysis']['total_probability']:.4f}")
    
    # Print off-target sites
    print("\nPotential Off-Target Sites:")
    for i, site in enumerate(results['sequence_analysis']['off_target_sites']):
        print(f"Site {i+1}: {site['sequence']} - Mismatches: {site['mismatches']} at positions {site['mismatch_positions']}")

