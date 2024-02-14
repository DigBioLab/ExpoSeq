import subprocess
import os

class WindowsInstall:    
    def start_wsl():
        filepath = os.path.join(os.path.dirname(__file__), "installation_windows.ps1")
        subprocess.run(["powershell.exe", "-File", filepath])
    
    def exit_wsl():
        subprocess.run(["wsl --shutdown"])
        
    def wsl():
        subprocess.run(["wsl"])
        
    @staticmethod
    def create_install_dir():
        subprocess.run(["wsl", "mkdir", "install_lcf"])
        
    @staticmethod
    def install_lcf():
        subprocess.run(["wsl", "fsutil", "file", "SetCaseSensitiveInfo", ".", "enable"])
        subprocess.run(["wsl", "wget", "https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh"])
        subprocess.run(["bash", "install_colabbatch_linux.sh"])
        

WindowsInstall.start_wsl()

    