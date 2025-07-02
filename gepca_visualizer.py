"""
Visualization Module for Gene Editing Precision Calculator Algorithm (GEPCA)

This module provides visualization tools for interpreting and presenting the results
of the Gene Editing Precision Calculator Algorithm.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

class GEPCAVisualizer:
    """
    Visualization class for the Gene Editing Precision Calculator Algorithm.
    """
    
    def __init__(self, theme: str = 'default'):
        """
        Initialize the visualizer.
        
        Args:
            theme: Visual theme for plots ('default', 'dark', or 'publication')
        """
        self.theme = theme
        self._set_theme()
    
    def _set_theme(self):
        """Set the visual theme for plots."""
        if self.theme == 'dark':
            plt.style.use('dark_background')
            self.colors = {
                'primary': '#00a8e8',
                'secondary': '#ff6b6b',
                'tertiary': '#98c379',
                'quaternary': '#c678dd',
                'background': '#282c34',
                'text': '#ffffff'
            }
        elif self.theme == 'publication':
            plt.style.use('seaborn-whitegrid')
            self.colors = {
                'primary': '#1f77b4',
                'secondary': '#ff7f0e',
                'tertiary': '#2ca02c',
                'quaternary': '#d62728',
                'background': '#ffffff',
                'text': '#000000'
            }
        else:  # default
            plt.style.use('seaborn')
            self.colors = {
                'primary': '#4c72b0',
                'secondary': '#dd8452',
                'tertiary': '#55a868',
                'quaternary': '#c44e52',
                'background': '#ffffff',
                'text': '#000000'
            }
    
    def plot_precision_score(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot the precision score and its components.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        precision_score = results['precision_score']
        components = precision_score['components']
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot precision score gauge
        self._plot_gauge(ax1, precision_score['precision_score'], 
                        "Precision Score", precision_score['confidence_interval'])
        
        # Plot components
        component_names = ['On-Target\nEfficiency', 'Off-Target\nProbability', 
                          'Structural\nContext', 'Edit\nCharacteristics']
        component_values = [
            components['on_target_efficiency'],
            components['off_target_probability'],
            components['structural_score'],
            components['edit_char_score']
        ]
        
        colors = [self.colors['primary'], self.colors['secondary'], 
                 self.colors['tertiary'], self.colors['quaternary']]
        
        ax2.bar(component_names, component_values, color=colors, alpha=0.7)
        ax2.set_ylim(0, 1)
        ax2.set_title('Precision Score Components')
        ax2.set_ylabel('Score')
        
        # Add horizontal line at 0.5 for reference
        ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def _plot_gauge(self, ax, value: float, title: str, confidence_interval: Tuple[float, float]):
        """
        Plot a gauge chart for the precision score.
        
        Args:
            ax: Matplotlib axis
            value: Value to plot (0-1)
            title: Title for the gauge
            confidence_interval: Tuple of (lower, upper) confidence bounds
        """
        # Define gauge properties
        gauge_min = 0
        gauge_max = 1
        
        # Create gauge
        theta = np.linspace(3*np.pi/4, 9*np.pi/4, 100)
        r = 0.5
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        
        # Plot gauge background
        ax.plot(x, y, color='gray', alpha=0.3, linewidth=10)
        
        # Calculate value position
        value_theta = 3*np.pi/4 + (9*np.pi/4 - 3*np.pi/4) * (value - gauge_min) / (gauge_max - gauge_min)
        value_x = r * np.cos(value_theta)
        value_y = r * np.sin(value_theta)
        
        # Plot confidence interval
        ci_lower_theta = 3*np.pi/4 + (9*np.pi/4 - 3*np.pi/4) * (confidence_interval[0] - gauge_min) / (gauge_max - gauge_min)
        ci_upper_theta = 3*np.pi/4 + (9*np.pi/4 - 3*np.pi/4) * (confidence_interval[1] - gauge_min) / (gauge_max - gauge_min)
        ci_theta = np.linspace(ci_lower_theta, ci_upper_theta, 50)
        ci_x = r * np.cos(ci_theta)
        ci_y = r * np.sin(ci_theta)
        ax.plot(ci_x, ci_y, color=self.colors['secondary'], alpha=0.5, linewidth=10)
        
        # Plot value
        ax.plot([0, value_x], [0, value_y], color=self.colors['primary'], linewidth=3)
        
        # Add a center circle
        circle = plt.Circle((0, 0), 0.05, color=self.colors['primary'])
        ax.add_artist(circle)
        
        # Add labels
        ax.text(-0.6, -0.6, '0.0', fontsize=12)
        ax.text(0.6, -0.6, '1.0', fontsize=12)
        
        # Add value text
        ax.text(0, -0.2, f"{value:.2f}", fontsize=24, ha='center', color=self.colors['primary'])
        
        # Add title
        ax.text(0, 0.7, title, fontsize=16, ha='center')
        
        # Add confidence interval text
        ax.text(0, -0.3, f"CI: [{confidence_interval[0]:.2f}, {confidence_interval[1]:.2f}]", 
               fontsize=10, ha='center', color=self.colors['secondary'])
        
        # Remove axis
        ax.set_xlim(-0.7, 0.7)
        ax.set_ylim(-0.7, 0.7)
        ax.axis('off')
    
    def plot_off_target_analysis(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot the off-target analysis results.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        off_target_analysis = results['off_target_analysis']
        sequence_analysis = results['sequence_analysis']
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot off-target probabilities
        site_ids = list(off_target_analysis['site_probabilities'].keys())
        probabilities = list(off_target_analysis['site_probabilities'].values())
        
        # Sort by probability
        sorted_indices = np.argsort(probabilities)[::-1]
        site_ids = [site_ids[i] for i in sorted_indices]
        probabilities = [probabilities[i] for i in sorted_indices]
        
        ax1.bar(site_ids, probabilities, color=self.colors['secondary'], alpha=0.7)
        ax1.set_title('Off-Target Site Probabilities')
        ax1.set_xlabel('Site ID')
        ax1.set_ylabel('Probability')
        ax1.set_ylim(0, max(probabilities) * 1.2)
        plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
        
        # Plot mismatch visualization
        off_target_sites = sequence_analysis['off_target_sites']
        guide_rna = sequence_analysis['target_site'].get('sequence', '')
        
        if guide_rna and off_target_sites:
            # Create a matrix for visualization
            num_sites = len(off_target_sites)
            guide_length = len(guide_rna)
            
            # Initialize matrix with zeros (matches)
            mismatch_matrix = np.zeros((num_sites, guide_length))
            
            # Fill in mismatches
            for i, site in enumerate(off_target_sites):
                for pos in site['mismatch_positions']:
                    if 0 <= pos < guide_length:
                        mismatch_matrix[i, pos] = 1
            
            # Plot heatmap
            sns.heatmap(mismatch_matrix, cmap=['#4CAF50', '#F44336'], 
                       cbar=False, ax=ax2, linewidths=0.5)
            
            # Add labels
            ax2.set_title('Mismatch Positions')
            ax2.set_xlabel('Position in Guide RNA')
            ax2.set_ylabel('Off-Target Site')
            ax2.set_yticks(np.arange(num_sites) + 0.5)
            ax2.set_yticklabels([f"Site {i+1}" for i in range(num_sites)])
            
            # Add legend
            legend_elements = [
                Patch(facecolor='#4CAF50', label='Match'),
                Patch(facecolor='#F44336', label='Mismatch')
            ]
            ax2.legend(handles=legend_elements, loc='upper right')
            
            # Add PAM region indicator
            ax2.axvline(x=guide_length - 3, color='black', linestyle='--', alpha=0.5)
            ax2.text(guide_length - 1.5, -0.5, 'PAM', ha='center')
            
            # Add seed region indicator
            ax2.axvline(x=guide_length - 12, color='black', linestyle=':', alpha=0.5)
            ax2.text(guide_length - 7.5, -0.5, 'Seed Region', ha='center')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_edit_patterns(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot the predicted edit patterns.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        editing_outcomes = results['editing_outcomes']
        edit_patterns = editing_outcomes['edit_patterns']
        
        # Create figure
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Plot pie chart of edit patterns
        labels = list(edit_patterns.keys())
        sizes = list(edit_patterns.values())
        
        # Custom colors for edit types
        edit_colors = {
            'insertion': self.colors['primary'],
            'deletion': self.colors['secondary'],
            'substitution': self.colors['tertiary'],
            'other': self.colors['quaternary']
        }
        
        colors = [edit_colors.get(label, self.colors['quaternary']) for label in labels]
        
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors,
              wedgeprops={'alpha': 0.7})
        ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
        ax.set_title('Predicted Edit Patterns')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_target_sequence_analysis(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot the target sequence analysis.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        sequence_analysis = results['sequence_analysis']
        target_site = sequence_analysis['target_site']
        
        if not target_site.get('found', False):
            print("Target site not found in sequence.")
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 3))
        
        # Get sequence and position
        sequence = target_site['sequence']
        position = target_site['position']
        
        # Plot sequence as colored blocks
        for i, base in enumerate(sequence):
            color = {
                'A': '#4CAF50',  # Green
                'T': '#F44336',  # Red
                'G': '#2196F3',  # Blue
                'C': '#FFC107'   # Yellow
            }.get(base, 'gray')
            
            ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, alpha=0.7))
            ax.text(i + 0.5, 0.5, base, ha='center', va='center', fontsize=12)
        
        # Add PAM site indicator (assuming PAM is after the guide RNA)
        if len(sequence) >= 23:  # If we have enough sequence to show PAM
            for i in range(3):
                if i + 20 < len(sequence):
                    ax.add_patch(plt.Rectangle((i + 20, 0), 1, 1, color='gray', alpha=0.3))
            ax.text(21.5, -0.3, 'PAM', ha='center')
        
        # Add seed region indicator
        ax.plot([8, 19], [-0.2, -0.2], 'k-', linewidth=2)
        ax.text(13.5, -0.5, 'Seed Region', ha='center')
        
        # Set axis limits
        ax.set_xlim(-0.5, len(sequence) + 0.5)
        ax.set_ylim(-1, 1.5)
        
        # Remove axis
        ax.axis('off')
        
        # Add title
        ax.set_title(f"Target Site Analysis (Position: {position}, Strand: {target_site['strand']})")
        
        # Add legend
        legend_elements = [
            Patch(facecolor='#4CAF50', label='A'),
            Patch(facecolor='#F44336', label='T'),
            Patch(facecolor='#2196F3', label='G'),
            Patch(facecolor='#FFC107', label='C')
        ]
        ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.1),
                 ncol=4)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_structural_context(self, results: Dict, save_path: Optional[str] = None):
        """
        Plot the structural context analysis.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        structural_context = results['structural_context']
        
        # Create figure
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Create radar chart
        categories = ['Chromatin\nAccessibility', 'DNA Shape\nScore', 'Structural\nScore']
        values = [
            structural_context['chromatin_accessibility'],
            structural_context['dna_shape_score'],
            structural_context['structural_score']
        ]
        
        # Add methylation status as text
        methylation_status = structural_context['methylation_status']
        
        # Number of variables
        N = len(categories)
        
        # What will be the angle of each axis in the plot (divide the plot / number of variables)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Values for the plot, plus close the loop
        values += values[:1]
        
        # Draw the plot
        ax.plot(angles, values, linewidth=2, linestyle='solid', color=self.colors['primary'])
        ax.fill(angles, values, color=self.colors['primary'], alpha=0.25)
        
        # Fix axis to go in the right order and start at 12 o'clock
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        
        # Draw axis lines for each angle and label
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        
        # Draw y-axis lines
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'])
        ax.set_ylim(0, 1)
        
        # Add title and methylation status
        ax.set_title('Structural Context Analysis')
        fig.text(0.5, 0.02, f"Methylation Status: {methylation_status}", ha='center')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def create_summary_dashboard(self, results: Dict, save_path: Optional[str] = None):
        """
        Create a comprehensive dashboard of all analysis results.
        
        Args:
            results: Results dictionary from the GEPCA analysis
            save_path: Optional path to save the figure
        """
        # Create figure with subplots
        fig = plt.figure(figsize=(15, 12))
        
        # Define grid layout
        gs = fig.add_gridspec(3, 3)
        
        # Precision score gauge (top left)
        ax1 = fig.add_subplot(gs[0, 0])
        precision_score = results['precision_score']
        self._plot_gauge(ax1, precision_score['precision_score'], 
                        "Precision Score", precision_score['confidence_interval'])
        
        # Precision components (top middle)
        ax2 = fig.add_subplot(gs[0, 1])
        components = precision_score['components']
        component_names = ['On-Target\nEff.', 'Off-Target\nProb.', 
                          'Structural\nContext', 'Edit\nChar.']
        component_values = [
            components['on_target_efficiency'],
            components['off_target_probability'],
            components['structural_score'],
            components['edit_char_score']
        ]
        colors = [self.colors['primary'], self.colors['secondary'], 
                 self.colors['tertiary'], self.colors['quaternary']]
        ax2.bar(component_names, component_values, color=colors, alpha=0.7)
        ax2.set_ylim(0, 1)
        ax2.set_title('Precision Components')
        
        # Edit patterns (top right)
        ax3 = fig.add_subplot(gs[0, 2])
        edit_patterns = results['editing_outcomes']['edit_patterns']
        labels = list(edit_patterns.keys())
        sizes = list(edit_patterns.values())
        edit_colors = {
            'insertion': self.colors['primary'],
            'deletion': self.colors['secondary'],
            'substitution': self.colors['tertiary'],
            'other': self.colors['quaternary']
        }
        colors = [edit_colors.get(label, self.colors['quaternary']) for label in labels]
        ax3.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors,
               wedgeprops={'alpha': 0.7})
        ax3.axis('equal')
        ax3.set_title('Edit Patterns')
        
        # Target sequence (middle row, span all columns)
        ax4 = fig.add_subplot(gs[1, :])
        sequence_analysis = results['sequence_analysis']
        target_site = sequence_analysis['target_site']
        
        if target_site.get('found', False):
            sequence = target_site['sequence']
            position = target_site['position']
            
            for i, base in enumerate(sequence):
                color = {
                    'A': '#4CAF50',
                    'T': '#F44336',
                    'G': '#2196F3',
                    'C': '#FFC107'
                }.get(base, 'gray')
                
                ax4.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, alpha=0.7))
                ax4.text(i + 0.5, 0.5, base, ha='center', va='center', fontsize=12)
            
            # Add PAM site indicator
            if len(sequence) >= 23:
                for i in range(3):
                    if i + 20 < len(sequence):
                        ax4.add_patch(plt.Rectangle((i + 20, 0), 1, 1, color='gray', alpha=0.3))
                ax4.text(21.5, -0.3, 'PAM', ha='center')
            
            # Add seed region indicator
            ax4.plot([8, 19], [-0.2, -0.2], 'k-', linewidth=2)
            ax4.text(13.5, -0.5, 'Seed Region', ha='center')
            
            ax4.set_xlim(-0.5, len(sequence) + 0.5)
            ax4.set_ylim(-1, 1.5)
            ax4.axis('off')
            ax4.set_title(f"Target Site (Position: {position}, Strand: {target_site['strand']})")
            
            # Add legend
            legend_elements = [
                Patch(facecolor='#4CAF50', label='A'),
                Patch(facecolor='#F44336', label='T'),
                Patch(facecolor='#2196F3', label='G'),
                Patch(facecolor='#FFC107', label='C')
            ]
            ax4.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.1),
                     ncol=4)
        else:
            ax4.text(0.5, 0.5, "Target site not found in sequence", 
                    ha='center', va='center', fontsize=14)
            ax4.axis('off')
        
        # Off-target probabilities (bottom left)
        ax5 = fig.add_subplot(gs[2, 0])
        off_target_analysis = results['off_target_analysis']
        site_ids = list(off_target_analysis['site_probabilities'].keys())
        probabilities = list(off_target_analysis['site_probabilities'].values())
        
        # Sort by probability
        sorted_indices = np.argsort(probabilities)[::-1]
        site_ids = [site_ids[i] for i in sorted_indices]
        probabilities = [probabilities[i] for i in sorted_indices]
        
        # Shorten site IDs for display
        short_ids = [f"Site {i+1}" for i in range(len(site_ids))]
        
        ax5.bar(short_ids, probabilities, color=self.colors['secondary'], alpha=0.7)
        ax5.set_title('Off-Target Probabilities')
        ax5.set_ylim(0, max(probabilities) * 1.2)
        plt.setp(ax5.get_xticklabels(), rotation=45, ha='right')
        
        # Structural context radar (bottom middle)
        ax6 = fig.add_subplot(gs[2, 1], polar=True)
        structural_context = results['structural_context']
        
        categories = ['Chromatin\nAccess.', 'DNA\nShape', 'Structural\nScore']
        values = [
            structural_context['chromatin_accessibility'],
            structural_context['dna_shape_score'],
            structural_context['structural_score']
        ]
        
        N = len(categories)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]
        values += values[:1]
        
        ax6.plot(angles, values, linewidth=2, linestyle='solid', color=self.colors['primary'])
        ax6.fill(angles, values, color=self.colors['primary'], alpha=0.25)
        
        ax6.set_theta_offset(np.pi / 2)
        ax6.set_theta_direction(-1)
        ax6.set_xticks(angles[:-1])
        ax6.set_xticklabels(categories)
        ax6.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax6.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=8)
        ax6.set_ylim(0, 1)
        ax6.set_title('Structural Context')
        
        # Mismatch visualization (bottom right)
        ax7 = fig.add_subplot(gs[2, 2])
        off_target_sites = sequence_analysis['off_target_sites']
        guide_rna = target_site.get('sequence', '')
        
        if guide_rna and off_target_sites:
            num_sites = len(off_target_sites)
            guide_length = len(guide_rna)
            
            mismatch_matrix = np.zeros((num_sites, guide_length))
            
            for i, site in enumerate(off_target_sites):
                for pos in site['mismatch_positions']:
                    if 0 <= pos < guide_length:
                        mismatch_matrix[i, pos] = 1
            
            sns.heatmap(mismatch_matrix, cmap=['#4CAF50', '#F44336'], 
                       cbar=False, ax=ax7, linewidths=0.5)
            
            ax7.set_title('Mismatch Positions')
            ax7.set_xlabel('Position in Guide RNA')
            ax7.set_ylabel('Off-Target Site')
            ax7.set_yticks(np.arange(num_sites) + 0.5)
            ax7.set_yticklabels([f"Site {i+1}" for i in range(num_sites)])
            
            # Add PAM and seed region indicators
            ax7.axvline(x=guide_length - 3, color='black', linestyle='--', alpha=0.5)
            ax7.text(guide_length - 1.5, -0.5, 'PAM', ha='center')
            ax7.axvline(x=guide_length - 12, color='black', linestyle=':', alpha=0.5)
            ax7.text(guide_length - 7.5, -0.5, 'Seed', ha='center')
        else:
            ax7.text(0.5, 0.5, "No off-target data available", 
                    ha='center', va='center', fontsize=12)
            ax7.axis('off')
        
        # Add title to the figure
        fig.suptitle('Gene Editing Precision Calculator - Analysis Dashboard', fontsize=16)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for the suptitle
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig


# Example usage
if __name__ == "__main__":
    # Import the calculator
    from gepca_implementation import GeneEditingPrecisionCalculator
    
    # Initialize the calculator
    calculator = GeneEditingPrecisionCalculator(reference_genome="hg38.fa")
    
    # Example target sequence and guide RNA
    target_sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    guide_rna = "ACGTACGTACGTACGTACGT"
    
    # Perform analysis
    results = calculator.analyze(target_sequence, guide_rna)
    
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
    
    print("Visualizations created and saved.")

