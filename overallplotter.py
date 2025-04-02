import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

# Bck end is non interactive
plt.switch_backend('Agg')

class PlottingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Plot Selection Tool")
        self.root.geometry("450x550")
        self.root.resizable(False, False)
        
        # For clean exit
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # Create style so it looks consistent
        style = ttk.Style()
        style.configure("TButton", font=("Arial", 10))
        style.configure("TLabel", font=("Arial", 10))
        style.configure("TCheckbutton", font=("Arial", 10))
        
        # Create a frame for the input
        input_frame = ttk.LabelFrame(root, text="Plot Selection")
        input_frame.pack(padx=10, pady=10, fill="both", expand=True)
        
        # Prefix entry
        ttk.Label(input_frame, text="File Prefix:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.prefix_var = tk.StringVar()
        self.prefix_entry = ttk.Entry(input_frame, textvariable=self.prefix_var, width=30)
        self.prefix_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        
        # Directory button
        self.dir_var = tk.StringVar(value=os.getcwd())
        dir_button = ttk.Button(input_frame, text="Browse", command=self.browse_directory)
        dir_button.grid(row=0, column=2, padx=5, pady=5)
        
        # Directory label
        ttk.Label(input_frame, text="Input Directory:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.dir_label = ttk.Label(input_frame, textvariable=self.dir_var, font=("Arial", 8))
        self.dir_label.grid(row=1, column=1, columnspan=2, sticky="w", padx=5, pady=5)
        
        # Expected frequency for spectral leakage
        ttk.Label(input_frame, text="Expected Frequency (Hz):").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.expected_freq_var = tk.StringVar(value="0.00")
        self.expected_freq_entry = ttk.Entry(input_frame, textvariable=self.expected_freq_var, width=10)
        self.expected_freq_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")
        
        # Output directory
        ttk.Label(input_frame, text="Output Directory:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.output_dir_var = tk.StringVar(value=os.getcwd())
        self.output_dir_label = ttk.Label(input_frame, textvariable=self.output_dir_var, font=("Arial", 8))
        self.output_dir_label.grid(row=3, column=1, sticky="w", padx=5, pady=5)
        output_dir_button = ttk.Button(input_frame, text="Browse", command=self.browse_output_directory)
        output_dir_button.grid(row=3, column=2, padx=5, pady=5)
        
        # Separator
        ttk.Separator(input_frame, orient="horizontal").grid(row=4, column=0, columnspan=3, sticky="ew", pady=10)
        
        # Checkboxes for plot selection
        ttk.Label(input_frame, text="Select Plots to Generate:").grid(row=5, column=0, columnspan=3, sticky="w", padx=5, pady=5)
        
        # Create checkboxes for each plot type
        self.plot_vars = {
            "time_domain": tk.BooleanVar(value=True),
            "dft_components": tk.BooleanVar(value=True),
            "dft_magnitude": tk.BooleanVar(value=True),
            "window_function": tk.BooleanVar(value=True),
            "spectral_leakage": tk.BooleanVar(value=True)
        }
        
        # Plot descriptions
        plot_descriptions = {
            "time_domain": "Time-Domain Signal",
            "dft_components": "DFT Real & Imaginary vs Sample Index",
            "dft_magnitude": "DFT Magnitude vs. Frequency",
            "window_function": "Window Function",
            "spectral_leakage": "Spectral Leakage Analysis"
        }
        
        # Add checkboxes to the frame
        for i, (key, description) in enumerate(plot_descriptions.items()):
            ttk.Checkbutton(input_frame, text=description, variable=self.plot_vars[key]).grid(
                row=i+6, column=0, columnspan=3, sticky="w", padx=20, pady=2
            )
        
        # Buttons frame
        button_frame = ttk.Frame(root)
        button_frame.pack(pady=10, padx=10, fill="x")
        
        # Generate button
        generate_button = ttk.Button(button_frame, text="Generate Plots", command=self.generate_plots)
        generate_button.pack(side="right", padx=5)
        
        # Cancel button
        cancel_button = ttk.Button(button_frame, text="Cancel", command=root.destroy)
        cancel_button.pack(side="right", padx=5)
        
        # Status label
        self.status_var = tk.StringVar()
        self.status_label = ttk.Label(root, textvariable=self.status_var, foreground="blue")
        self.status_label.pack(pady=5)
    
    def browse_directory(self):
        """Allow user to select input directory"""
        directory = filedialog.askdirectory(initialdir=self.dir_var.get())
        if directory:  # If user didnt cancel
            self.dir_var.set(directory)
    
    def browse_output_directory(self):
        """Allow user to select output directory"""
        directory = filedialog.askdirectory(initialdir=self.output_dir_var.get())
        if directory:  # If user didnt cancel
            self.output_dir_var.set(directory)
    
    def on_closing(self):
        """Handle window close event"""
        plt.close('all')  # Close any open matplotlib figures
        self.root.destroy()
            
    def generate_plots(self):
        prefix = self.prefix_var.get().strip()
        
        if not prefix:
            messagebox.showerror("Error", "Please enter a file prefix")
            return
        
        # Get input and output directories
        input_dir = self.dir_var.get()
        output_dir = self.output_dir_var.get()
        
        # Make sure directories exist
        if not os.path.isdir(input_dir):
            messagebox.showerror("Error", f"Input directory does not exist: {input_dir}")
            return
            
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                messagebox.showerror("Error", f"Could not create output directory: {str(e)}")
                return
        
        # Get the expected frequency
        try:
            expected_freq = float(self.expected_freq_var.get())
        except ValueError:
            messagebox.showerror("Error", "Expected frequency must be a number")
            return
        
        # Get the selected plots
        selected_plots = {k: v.get() for k, v in self.plot_vars.items()}
        
        if not any(selected_plots.values()):
            messagebox.showerror("Error", "Please select at least one plot type")
            return
        
        try:
            # Generate the selected plots
            self.status_var.set("Generating plots...")
            self.root.update()
            
            # Time domain signal
            if selected_plots["time_domain"]:
                self.plot_time_domain(prefix)
            
            # DFT components
            if selected_plots["dft_components"]:
                self.plot_dft_components(prefix)
            
            # DFT magnitude
            if selected_plots["dft_magnitude"]:
                self.plot_dft_magnitude(prefix)
            
            # Window function
            if selected_plots["window_function"]:
                self.plot_window_function(prefix)
            
            # Spectral leakage
            if selected_plots["spectral_leakage"]:
                self.plot_spectral_leakage(prefix, expected_freq)
            
            self.status_var.set(f"All selected plots generated successfully!")
        except Exception as e:
            self.status_var.set(f"Error: {str(e)}")
            messagebox.showerror("Error", f"An error occurred: {str(e)}")
    
    def check_columns(self, df, required_columns, file_type):
        """Check if dataframe has required columns"""
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing columns in {file_type} file: {', '.join(missing_columns)}")
        return True
            
    def plot_time_domain(self, prefix):
        try:
            input_file = os.path.join(self.dir_var.get(), f"{prefix}_signal_raw.csv")
            output_file = os.path.join(self.output_dir_var.get(), f"{prefix}_time_domain.png")
            
            # Check if file exists
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"Could not find time domain data file: {input_file}")
                
            # Read data
            time_data = pd.read_csv(input_file, encoding='utf-8')
            
            # Check for required columns
            self.check_columns(time_data, ['time', 'real'], 'time domain')
            
            # Create plot
            fig = plt.figure(figsize=(8, 6), dpi=100)
            try:
                plt.plot(time_data['time'], time_data['real'], 'b-', linewidth=1.5)
                plt.title('Time-Domain Signal', fontsize=12)
                plt.xlabel('Time', fontsize=10)
                plt.ylabel('Amplitude', fontsize=10)
                plt.grid(False)  # Explicitly disable grid
                plt.tight_layout()
                plt.savefig(output_file, dpi=100)
            finally:
                plt.close(fig)  # Ensure figure is closed even if error occurs
                
            self.status_var.set(f"Generated: {os.path.basename(output_file)}")
            self.root.update()
            
        except Exception as e:
            raise Exception(f"Error in time domain plot: {str(e)}")
    
    def plot_dft_components(self, prefix):
        try:
            input_file = os.path.join(self.dir_var.get(), f"{prefix}_dft.csv")
            output_file = os.path.join(self.output_dir_var.get(), f"{prefix}_dft_components.png")
            
            # Check if file exists
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"Could not find DFT components data file: {input_file}")
                
            # Read data
            dft_data = pd.read_csv(input_file, encoding='utf-8')
            
            # Check for required columns
            self.check_columns(dft_data, ['index', 'real', 'imaginary'], 'DFT components')
            
            # Create plot
            fig = plt.figure(figsize=(8, 6), dpi=100)
            try:
                plt.plot(dft_data['index'], dft_data['real'], 'b-', label='Real', linewidth=1.5)
                plt.plot(dft_data['index'], dft_data['imaginary'], 'r-', label='Imaginary', linewidth=1.5)
                plt.title('DFT vs. Sample Index', fontsize=12)
                plt.xlabel('Sample Index', fontsize=10)
                plt.ylabel('DFT Value', fontsize=10)
                plt.legend(loc='best')
                plt.grid(False)  
                plt.tight_layout()
                plt.savefig(output_file, dpi=100)
            finally:
                plt.close(fig)
                
            self.status_var.set(f"Generated: {os.path.basename(output_file)}")
            self.root.update()
            
        except Exception as e:
            raise Exception(f"Error in DFT components plot: {str(e)}")
    
    def plot_dft_magnitude(self, prefix):
        try:
            input_file = os.path.join(self.dir_var.get(), f"{prefix}_magnitude.csv")
            output_file = os.path.join(self.output_dir_var.get(), f"{prefix}_dft_magnitude.png")
            
            # Check if file exists
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"Could not find magnitude data file: {input_file}")
                
            # Read data
            mag_data = pd.read_csv(input_file, encoding='utf-8')
            
            # Check for required columns
            self.check_columns(mag_data, ['frequency', 'magnitude'], 'DFT magnitude')
            
            # Create plot
            fig = plt.figure(figsize=(8, 6), dpi=100)
            try:
                # Add a small value to avoid log(0) issues
                epsilon = 1e-15
                mag_data['magnitude'] = mag_data['magnitude'] + epsilon
                plt.plot(mag_data['frequency'], mag_data['magnitude'], 'b-', linewidth=1.5)
                plt.title('DFT Magnitude vs. Frequency', fontsize=12)
                plt.xlabel('Frequency (Hz)', fontsize=10)
                plt.ylabel('Magnitude', fontsize=10)
                # Keep grid only for DFT magnitude as in the original code
                plt.grid(True, which="both", ls="-", alpha=0.3)
                plt.tight_layout()
                plt.savefig(output_file, dpi=100)
            finally:
                plt.close(fig)
                
            self.status_var.set(f"Generated: {os.path.basename(output_file)}")
            self.root.update()
            
        except Exception as e:
            raise Exception(f"Error in DFT magnitude plot: {str(e)}")
    
    def plot_window_function(self, prefix):
        try:
            input_file = os.path.join(self.dir_var.get(), f"{prefix}_window.csv")
            output_file = os.path.join(self.output_dir_var.get(), f"{prefix}_window_function.png")
            
            # Check if file exists
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"Could not find window function data file: {input_file}")
                
            # Read data
            window = pd.read_csv(input_file, encoding='utf-8')
            
            # Check for required columns
            self.check_columns(window, ['x', 'window_value'], 'window function')
            
            # Create plot
            fig = plt.figure(figsize=(10, 6), dpi=100)
            try:
                plt.plot(window['x'], window['window_value'], 'b-', linewidth=1.5)
                plt.title('Window Function', fontsize=12)
                plt.xlabel('Normalised Position', fontsize=10)
                plt.ylabel('Weight', fontsize=10)
                plt.grid(False)  # Explicitly disable grid
                plt.tight_layout()
                plt.savefig(output_file, dpi=100)
            finally:
                plt.close(fig)
                
            self.status_var.set(f"Generated: {os.path.basename(output_file)}")
            self.root.update()
            
        except Exception as e:
            raise Exception(f"Error in window function plot: {str(e)}")
    
    def plot_spectral_leakage(self, prefix, expected_freq=0.00):
        try:
            input_file = os.path.join(self.dir_var.get(), f"{prefix}_magnitude.csv")
            output_file = os.path.join(self.output_dir_var.get(), f"{prefix}_spectral_leakage.png")
            
            # Check if file exists
            if not os.path.isfile(input_file):
                raise FileNotFoundError(f"Could not find magnitude data file: {input_file}")
                
            # Read data
            data = pd.read_csv(input_file, encoding='utf-8')
            
            # Check for required columns
            self.check_columns(data, ['frequency', 'magnitude'], 'spectral leakage')
            
            # Find the actual peak
            peak_idx = data['magnitude'].idxmax()
            peak_freq = data.loc[peak_idx, 'frequency']
            peak_mag = data.loc[peak_idx, 'magnitude']
            
            # Create figure with two subplots
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), dpi=100)
            
            try:
                # Plot 1: Full spectrum view
                full_view_data = data[(data['frequency'] >= -5) & (data['frequency'] <= 5)]
                if full_view_data.empty:
                    full_view_data = data  # Fallback if no data in specified range
                
                ax1.plot(full_view_data['frequency'], full_view_data['magnitude'], 'b-', linewidth=1.5)
                ax1.set_title('Full Spectrum (-5-5 Hz)', fontsize=12)
                ax1.set_xlabel('Frequency (Hz)', fontsize=10)
                ax1.set_ylabel('Magnitude', fontsize=10)
                ax1.plot(peak_freq, peak_mag, 'ro', markersize=8,
                       label=f'Actual Peak ({peak_freq:.4f} Hz)')
                
                # Add expected frequency line if it's within range
                if -5 <= expected_freq <= 5:
                    ax1.axvline(x=expected_freq, color='g', linestyle='--', linewidth=1.5,
                              label=f'Expected ({expected_freq:.2f} Hz)')
                
                ax1.legend(loc='best')
                ax1.grid(False)  # Explicitly disable grid
                
                # Plot 2: Zoomed view around the peak
                zoom_min = max(min(data['frequency']), peak_freq - 2.0)
                zoom_max = min(max(data['frequency']), peak_freq + 2.0)
                
                # Get data in the zoom range
                zoom_data = data[(data['frequency'] >= zoom_min) & (data['frequency'] <= zoom_max)]
                if zoom_data.empty:
                    zoom_data = data  # Fallback if no data in specified range
                    
                ax2.plot(zoom_data['frequency'], zoom_data['magnitude'], 'b-', linewidth=2)
                ax2.set_title(f'Zoomed View ({zoom_min:.1f}-{zoom_max:.1f} Hz)', fontsize=12)
                ax2.set_xlabel('Frequency (Hz)', fontsize=10)
                ax2.set_ylabel('Magnitude', fontsize=10)
                ax2.plot(peak_freq, peak_mag, 'ro', markersize=8,
                       label=f'Actual Peak ({peak_freq:.4f} Hz)')
                
                # Add expected frequency line if it's within zoom range
                if zoom_min <= expected_freq <= zoom_max:
                    ax2.axvline(x=expected_freq, color='g', linestyle='--', linewidth=1.5,
                              label=f'Expected ({expected_freq:.2f} Hz)')
                
                # Annotate significant peaks in zoomed view
                significant_peaks = zoom_data[zoom_data['magnitude'] > peak_mag * 0.1]
                for _, row in significant_peaks.iterrows():
                    percentage = (row['magnitude'] / peak_mag) * 100
                    ax2.annotate(f"{percentage:.1f}%",
                               xy=(row['frequency'], row['magnitude']),
                               xytext=(5, 5),
                               textcoords='offset points',
                               fontsize=8)
                
                ax2.legend(loc='best')
                ax2.grid(False)  # Explicitly disable grid
                
                # Adjust layout and save
                plt.tight_layout()
                plt.savefig(output_file, dpi=100)
            finally:
                plt.close(fig)
                
            self.status_var.set(f"Generated: {os.path.basename(output_file)}")
            self.root.update()
            
        except Exception as e:
            raise Exception(f"Error in spectral leakage plot: {str(e)}")

# Can be used as a module incase of reuse in future
if __name__ == "__main__":
    root = tk.Tk()
    app = PlottingApp(root)
    root.mainloop()