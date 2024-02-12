# Check if running on Windows
if ($env:OS -ne "Windows_NT") {
    Write-Host "This script is intended to run on Windows."
    exit
}

# Enable WSL2
wsl --set-default-version 2

# Install a Linux Distribution (e.g., Ubuntu)
$distro = "Ubuntu"
$distroAppx = Get-AppxPackage | Where-Object { $_.Name -like "*$distro*" }
if ($distroAppx -eq $null) {
    Write-Host "Installing $distro..."
    Start-Process -FilePath "https://aka.ms/wsl$([System.Environment]::OSVersion.Version.Major)" -Wait
}

# Launch PowerShell as Administrator for symlink creation
#Start-Process -FilePath "powershell.exe" -ArgumentList "Start-Process -FilePath 'fsutil' -ArgumentList 'file SetCaseSensitiveInfo path\to\localcolabfold\installation enable' -Verb RunAs" -Wait

# Launch the installed Linux distribution
wsl

# Get the current directory
$currentDirectory = (Get-Location).Path

# Define the directory name
$directoryName = "localcolabfold"

# Combine the current directory path with the directory name
$directoryPath = Join-Path -Path $currentDirectory -ChildPath $directoryName

# Create the directory
New-Item -ItemType Directory -Path $directoryPath -Force

# Change directory to the newly created directory
cd $directoryPath

# Download the installation script using Invoke-WebRequest
Invoke-WebRequest -Uri "https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh" -OutFile "install_colabbatch_linux.sh"


wsl -e sh install_colabbatch_linux.sh

# Setup Environment Variables in Linux
$bashrcPath = "~/.bashrc"
Add-Content -Path $bashrcPath -Value "export TF_FORCE_UNIFIED_MEMORY='1'"
Add-Content -Path $bashrcPath -Value "export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'"
Add-Content -Path $bashrcPath -Value "export XLA_PYTHON_CLIENT_ALLOCATOR='platform'"
Add-Content -Path $bashrcPath -Value "export TF_FORCE_GPU_ALLOW_GROWTH='true'"



wsl --shutdown
Write-Host "Setup complete. Please restart your WSL terminal to apply the changes."
