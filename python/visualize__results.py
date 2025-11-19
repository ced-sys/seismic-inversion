"""
Visualize C++ Inversion Results
Reads the output from the C++ tool and creates publication-quality plots
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def visualize_inversion_results():
    """Create comprehensive visualization of inversion results"""
    
    print("=" * 60)
    print("INVERSION RESULTS VISUALIZATION")
    print("=" * 60)
    
    # Read results
    results_file = "../output/inversion_results.csv"
    
    if not os.path.exists(results_file):
        print(f"✗ Error: Results file not found: {results_file}")
        print("  Please run the C++ inversion tool first!")
        return
    
    print(f"\n[1/2] Loading results from {results_file}...")
    df = pd.read_csv(results_file)
    print(f"✓ Loaded {len(df)} samples")
    
    # Create figure with multiple subplots
    print("\n[2/2] Creating visualization...")
    fig = plt.figure(figsize=(16, 10))
    
    # Define time axis (in seconds)
    time = df['sample'].values / 20.0  # Assuming ~20 Hz sampling
    
    # Plot 1: Original vs Synthetic Trace
    ax1 = plt.subplot(4, 1, 1)
    ax1.plot(time, df['original_trace'], 'b-', linewidth=1, label='Original Trace', alpha=0.7)
    ax1.plot(time, df['synthetic_trace'], 'r--', linewidth=1.5, label='Synthetic Trace')
    ax1.set_ylabel('Amplitude')
    ax1.set_title('Seismic Inversion Results - Real IRIS Data', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([time[0], time[-1]])
    
    # Plot 2: Residual (Error)
    ax2 = plt.subplot(4, 1, 2)
    residual = df['original_trace'] - df['synthetic_trace']
    ax2.plot(time, residual, 'g-', linewidth=0.8, alpha=0.7)
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax2.set_ylabel('Residual')
    ax2.set_title('Inversion Residual (Original - Synthetic)')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([time[0], time[-1]])
    
    # Add RMS error text
    rms_error = np.sqrt(np.mean(residual**2))
    ax2.text(0.02, 0.95, f'RMS Error: {rms_error:.4f}', 
             transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', 
             facecolor='wheat', alpha=0.5))
    
    # Plot 3: Reflectivity Series
    ax3 = plt.subplot(4, 1, 3)
    # Use stem plot for reflectivity to show spikes clearly
    significant_spikes = np.abs(df['reflectivity']) > 0.01
    spike_times = time[significant_spikes]
    spike_values = df['reflectivity'][significant_spikes]
    
    ax3.stem(spike_times, spike_values, linefmt='b-', markerfmt='bo', 
             basefmt='k-', label='Reflectivity Spikes')
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=0.8)
    ax3.set_ylabel('Reflectivity')
    ax3.set_title(f'Reflectivity Series ({len(spike_times)} significant spikes)')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([time[0], time[-1]])
    ax3.legend(loc='upper right')
    
    # Plot 4: Acoustic Impedance
    ax4 = plt.subplot(4, 1, 4)
    ax4.plot(time, df['impedance'], 'r-', linewidth=2, label='Acoustic Impedance')
    ax4.fill_between(time, df['impedance'], alpha=0.3, color='red')
    ax4.set_xlabel('Time (seconds)', fontsize=11)
    ax4.set_ylabel('Relative Impedance')
    ax4.set_title('Acoustic Impedance Profile')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([time[0], time[-1]])
    ax4.legend(loc='upper right')
    
    # Add statistics box
    stats_text = f'Min: {df["impedance"].min():.3f}\n'
    stats_text += f'Max: {df["impedance"].max():.3f}\n'
    stats_text += f'Mean: {df["impedance"].mean():.3f}'
    ax4.text(0.02, 0.95, stats_text, transform=ax4.transAxes, 
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    plt.tight_layout()
    
    # Save figure
    output_file = "../output/inversion_visualization.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    
    # Create a second figure with zoomed views
    fig2, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Find interesting regions (high activity)
    window_size = 1000
    activity = np.convolve(np.abs(df['reflectivity']), 
                          np.ones(window_size)/window_size, mode='same')
    interesting_idx = np.argmax(activity)
    
    zoom_start = max(0, interesting_idx - window_size//2)
    zoom_end = min(len(df), interesting_idx + window_size//2)
    
    zoom_time = time[zoom_start:zoom_end]
    
    # Zoomed Plot 1: Traces
    axes[0, 0].plot(zoom_time, df['original_trace'][zoom_start:zoom_end], 
                    'b-', linewidth=1.5, label='Original', alpha=0.7)
    axes[0, 0].plot(zoom_time, df['synthetic_trace'][zoom_start:zoom_end], 
                    'r--', linewidth=2, label='Synthetic')
    axes[0, 0].set_title('Zoomed: Trace Comparison')
    axes[0, 0].set_ylabel('Amplitude')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Zoomed Plot 2: Reflectivity
    zoom_spikes = df['reflectivity'][zoom_start:zoom_end]
    significant = np.abs(zoom_spikes) > 0.01
    axes[0, 1].stem(zoom_time[significant], zoom_spikes[significant], 
                    linefmt='b-', markerfmt='bo', basefmt='k-')
    axes[0, 1].set_title('Zoomed: Reflectivity Detail')
    axes[0, 1].set_ylabel('Reflectivity')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Zoomed Plot 3: Impedance
    axes[1, 0].plot(zoom_time, df['impedance'][zoom_start:zoom_end], 
                    'r-', linewidth=2)
    axes[1, 0].fill_between(zoom_time, 
                           df['impedance'][zoom_start:zoom_end], 
                           alpha=0.3, color='red')
    axes[1, 0].set_title('Zoomed: Impedance Profile')
    axes[1, 0].set_xlabel('Time (seconds)')
    axes[1, 0].set_ylabel('Impedance')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Correlation Analysis
    correlation = np.correlate(df['original_trace'], df['synthetic_trace'], 
                              mode='same')
    correlation = correlation / np.max(np.abs(correlation))
    lag_time = (np.arange(len(correlation)) - len(correlation)//2) / 20.0
    axes[1, 1].plot(lag_time, correlation, 'purple', linewidth=2)
    axes[1, 1].axvline(x=0, color='k', linestyle='--', linewidth=1)
    axes[1, 1].set_title('Cross-Correlation (Fit Quality)')
    axes[1, 1].set_xlabel('Lag (seconds)')
    axes[1, 1].set_ylabel('Normalized Correlation')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_xlim([-50, 50])
    
    # Add correlation coefficient
    corr_coef = np.corrcoef(df['original_trace'], df['synthetic_trace'])[0, 1]
    axes[1, 1].text(0.02, 0.95, f'Correlation: {corr_coef:.4f}', 
                   transform=axes[1, 1].transAxes, fontsize=10,
                   verticalalignment='top', 
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    plt.tight_layout()
    output_file2 = "../output/inversion_details.png"
    plt.savefig(output_file2, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file2}")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Total samples:           {len(df)}")
    print(f"Duration:                {time[-1]:.2f} seconds")
    print(f"RMS Error:               {rms_error:.6f}")
    print(f"Correlation Coefficient: {corr_coef:.6f}")
    print(f"Number of spikes:        {len(spike_times)}")
    print(f"Impedance range:         {df['impedance'].min():.3f} to {df['impedance'].max():.3f}")
    print(f"Impedance contrast:      {(df['impedance'].max() - df['impedance'].min()) / df['impedance'].min() * 100:.2f}%")
    print("=" * 60)
    
    print("\n✓ Visualization complete!")
    print("\nGenerated files:")
    print(f"  1. {output_file}")
    print(f"  2. {output_file2}")
    print("\nYou can now include these figures in your report!\n")

if __name__ == "__main__":
    try:
        visualize_inversion_results()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
