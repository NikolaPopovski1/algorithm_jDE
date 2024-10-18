import subprocess
import re
import numpy as np

# Lists to store energy, runtime, nfes, and speed
energies = []
runtimes = []
nfes_list = []
speeds = []

# Run the program for 50 seeds
for seed in range(1, 51):
    result = subprocess.run([f'.\\a.exe', "ABBBBBBABBBAB", "-seed", str(seed), "-target", "-5.6104", 
                             "-nfesLmt", "1000000", "-Np", "300", "-runtimeLmt", "3600"], capture_output=True, text=True)
    
    # Extract relevant values using regex
    time_taken = float(re.search(r"Time taken: (\d+)", result.stdout).group(1))
    best_energy = float(re.search(r"Best energy: ([\d\.\-]+)", result.stdout).group(1))
    nfes_match = re.search(r"Speed.*nfes / duration = speed\): (\d+)", result.stdout)
    if nfes_match:
        nfes = int(nfes_match.group(1))
    else:
        print(f"Could not find nfes in output for seed {seed}")
    speed = float(re.search(r"Speed.*= ([\d\.\-]+)", result.stdout).group(1))

    # Append values to lists
    energies.append(best_energy)
    runtimes.append(time_taken)
    nfes_list.append(nfes)
    speeds.append(speed)

# Calculate statistics
best_energy = min(energies)
worst_energy = max(energies)
average_energy = np.mean(energies)
std_dev_energy = np.std(energies)
average_runtime = np.mean(runtimes)
average_nfes = np.mean(nfes_list)
average_speed = np.mean(speeds)

# Output results
print(f"Best energy: {best_energy}")
print(f"Worst energy: {worst_energy}")
print(f"Average energy: {average_energy}")
print(f"Standard deviation of energy: {std_dev_energy}")
print(f"Average runtime: {average_runtime} ms")
print(f"Average nfes: {average_nfes}")
print(f"Average speed: {average_speed}")