import os
import subprocess
import tarfile
import requests

class NetSolP:
    
    def launch():
        save_dir =  os.path.join(os.path.dirname(__file__), "NetSolP-1.0", "model_weights")
        os.makedirs(save_dir, exists_ok = True)
        response = requests.get(r"https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz")
        if response.status_code == 200:
            # Open the file and write the content in binary mode
            with open(os.path.join(), "wb") as f:
                f.write(response.content)
            print("File downloaded successfully and saved as:", filename)

            # Extract the contents of the downloaded file
            extract_dir = os.path.join(save_directory, "extracted")
            os.makedirs(extract_dir, exist_ok=True)
            with tarfile.open(filename, "r:gz") as tar:
                tar.extractall(extract_dir)

    
        