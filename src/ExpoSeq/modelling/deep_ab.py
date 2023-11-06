from .set_env import VE
import subprocess
import os

class DeepAB:
    def __init__(self):
        self.git_available()
        Ve_Control = VE(name_ve = "deepab")
        self.clone_repo()
        Ve_Control.switch_to_ve()
        self.install_requirements()
        
        
    def clone_repo():
        if not os.path.isdir("DeepAb"):
            while user_input not in ["X", "x", ""]:
                user_input = input("You will now clone the repository from DeepAb.\nPlease make sure that you have read the corresponding licencing on their repository.\nPress enter to continue and press x to abort the process.")
                if user_input not in ["X", "x", ""]:
                    print("Please enter a correct value.")
            
            if user_input == "": 
                subprocess.run(["git",
                                "clone", 
                                "https://github.com/RosettaCommons/DeepAb.git"])
            else:
                pass
            
    def install_requirements():
        if os.path.isfile(os.path.join(os.getcwd(), "DeepAb", "requirements.txt")):
            subprocess.run(["pip",
                            "install", 
                            "-r", 
                            "DeepAb/requirements.txt"])
            
    

