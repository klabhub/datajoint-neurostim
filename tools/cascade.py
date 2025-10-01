"""
Wrapper function for running CASCADE spike inference from Calcium signals stored in MATLAB files.

@author: Bart Krekelberg

"""

import os, sys, glob
import numpy as np
import scipy.io as sio
import ruamel.yaml as yaml
import argparse
yaml = yaml.YAML(typ='rt')

def run_cascade_inference(dffFile, model_name, 
                         model_folder="Pretrained_models", 
                         threshold=0, 
                         padding=np.nan, 
                         trace_noise_levels=None, 
                         verbosity=1,
                         cascade_folder=None):
    """
    Run CASCADE spike inference on dF/F data from MATLAB files.
    
    Parameters
    ----------
    dffFile : str
        Path to the file containing the input MATLAB file with dF/F data.
        (A .mat file with a variable 'dFF' of shape (timepoints, neurons))
    cascade_folder : str, optional
        Path to CASCADE installation folder (default: current working directory)
    model_name : str
        Name of the CASCADE model to use (e.g., 'Global_EXC_15Hz_smoothing100ms')
    model_folder : str, optional
        Path to folder containing CASCADE models (default: "Pretrained_models")
    threshold : int or bool, optional
        Thresholding option for spike probabilities (default: 0)
        0: Set negative values to 0
        1 or True: Apply spike threshold with dilated mask
        False: No thresholding
    padding : float, optional
        Value for padding where predictions cannot be made (default: np.nan)
    trace_noise_levels : array-like, optional
        Noise levels of traces. If None, calculated automatically (default: None)
    verbosity : int, optional
        Verbosity level (0 or 1, default: 1)
            
    Returns
    -------
    output_path : str
        Full path to saved output file
        
    Raises
    ------
    FileNotFoundError
        If input MATLAB file is not found
    KeyError
        If specified variable is not found in MATLAB file
    """
    
    # Setup CASCADE path - default to current working directory
    if cascade_folder is None:
        cascade_folder = os.getcwd()
    
    original_dir = os.getcwd()
    sys.path.append(cascade_folder)
    os.chdir(cascade_folder)
    
    try:
        from cascade2p import cascade
        from cascade2p.utils_discrete_spikes import infer_discrete_spikes
      
        input_variable = 'dFF'
        output_filename = 'results.mat'
    
        
        # Load dF/F data from MATLAB file
        if not os.path.exists(dffFile):
            raise FileNotFoundError(f"Input file not found: {dffFile}")
            
        if verbosity >= 1:
            print(f"Loading dF/F data from: {dffFile}")
            
        mat_data = sio.loadmat(dffFile)
        
        if input_variable not in mat_data:
            available_vars = [k for k in mat_data.keys() if not k.startswith('__')]
            raise KeyError(f"Variable '{input_variable}' not found in MATLAB file. "
                          f"Available variables: {available_vars}")
        
        dFF = mat_data[input_variable].T  # Transpose to shape (neurons, timepoints) for cascade
        
        if verbosity >= 1:
            print(f"Loaded dF/F data shape: {dFF.shape}")
            print(f"Using CASCADE model: {model_name}")
        
        # Check if model already exists before downloading
        model_path = os.path.join(model_folder, model_name)
        cfg_file = os.path.join(model_path, "config.yaml")
        
        if os.path.isfile(cfg_file):
            if verbosity >= 1:
                print(f"Model {model_name} already exists at {model_path}")
        else:
            if verbosity >= 1:
                print(f"Model {model_name} not found. Downloading...")
            cascade.download_model(model_name, model_folder=model_folder, verbose=verbosity)
        
        # Run spike inference
        spike_prob = cascade.predict(
            model_name=model_name,
            traces=dFF,
            model_folder=model_folder,
            threshold=threshold,
            padding=padding,
            trace_noise_levels=trace_noise_levels,
            verbosity=verbosity
        )
        
        discrete, spikeSample = infer_discrete_spikes(
            spike_prob, model_name, 
            model_folder=model_folder, 
            verbosity=verbosity)
        
        # Convert spikeSample to cell array format for MATLAB
        # Each cell contains spike samples as a column vector
        # Handle both matrix and list formats from infer_discrete_spikes
        spikeSample_cell = []
        if isinstance(spikeSample, np.ndarray):
            # Convert matrix to list of arrays (one per neuron)
            spikeSample = [spikeSample[i, :] for i in range(spikeSample.shape[0])]
        for neuron_spikes in spikeSample:
            if len(neuron_spikes) > 0:
                spikeSample_cell.append(neuron_spikes.reshape(-1, 1))  # Column vector
            else:
                spikeSample_cell.append(np.array([]).reshape(0, 1))   # Empty column vector
        spikeSample = spikeSample_cell
       
        
        # Save results
        # Create output path by replacing extension with _cascade.mat
        base_name = os.path.splitext(dffFile)[0]
        output_path = base_name + '_cascade.mat'
        sio.savemat(output_path, {'pSpike': spike_prob.T , 'sSpike': spikeSample})

        if verbosity >= 1:
            print(f"Results saved to: {output_path}")            

        return output_path

    finally:
        # Return to original directory
        os.chdir(original_dir)


def main():
    """
    Command line interface for CASCADE spike inference.
    """
    parser = argparse.ArgumentParser(description='Run CASCADE spike inference on calcium imaging data')
    
    # Required arguments
    parser.add_argument('dff_file', type=str, 
                       help='Path to MATLAB file containing dF/F data (variable "dFF")')
    parser.add_argument('model_name', type=str,
                       help='Name of CASCADE model (e.g., "Global_EXC_15Hz_smoothing200ms")')
    
    # Optional arguments
    parser.add_argument('--cascade_folder', type=str, 
                       default='C:/Users/bartk/OneDrive - Rutgers University/Documents/common/github/Cascade',
                       help='Path to CASCADE installation folder')
    parser.add_argument('--model_folder', type=str, default='Pretrained_models',
                       help='Path to folder containing CASCADE models (default: "Pretrained_models")')
    parser.add_argument('--threshold', type=int, default=0, choices=[0, 1],
                       help='Thresholding option: 0=negative to zero, 1=spike threshold (default: 0)')
    parser.add_argument('--verbosity', type=int, default=1, choices=[0, 1],
                       help='Verbosity level: 0=quiet, 1=verbose (default: 1)')
    
    args = parser.parse_args()
    
    # Run CASCADE inference
    output_path = run_cascade_inference(
        dffFile=args.dff_file,
        model_name=args.model_name,
        cascade_folder=args.cascade_folder,
        model_folder=args.model_folder,
        threshold=args.threshold,
        verbosity=args.verbosity
    )
    
    if args.verbosity >= 1:
        print("CASCADE inference complete!")
        print(f"Results saved to: {output_path}")


if __name__ == "__main__":
    main()
