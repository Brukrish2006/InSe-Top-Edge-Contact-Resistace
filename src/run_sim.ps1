# Set the DEVSIM Math Libraries path to Intel MKL
$env:DEVSIM_MATH_LIBS="C:\Python314\Library\bin\mkl_rt.3.dll"

# Make sure we add the python_packages from DEVSIM to PYTHONPATH just in case
$env:PYTHONPATH="C:\Users\harsh\devsim\python_packages"

Write-Host "DEVSIM environment is set."
Write-Host "Running InSe simulations..."

# Run the python script
python run_inse_simulations.py

Write-Host "Simulations complete. Results are saved in summary_results.csv."
