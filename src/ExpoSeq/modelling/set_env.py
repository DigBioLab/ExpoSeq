import subprocess
import os

class VE:
    
    def __init__(self, name_ve):
        self.name_ve
        self.create_ve()
        self.switch_to_ve()
        
    
    def check_path(self):
        wd = os.getcwd()
        ve_existence = os.path.isdir(os.path.join(wd, self.name_ve))
        return ve_existence
    
    def check_venv_correct(self):
        wd = os.getcwd()
        exists = os.path.isfile(os.path.join(wd, self.name_ve, "pyvenv.cfg"))
        return exists
    
    def create_ve(self):
        if self.check_path() == False:
            subprocess.run(["python",
                        "-m",
                        "venv",
                        self.name_ve])
        elif self.check_venv_correct() == False:
            subprocess.run(["python",
                            "-m",
                            "venv",
                            self.name_ve + "_ve"])        
    def switch_to_ve(self):
        subprocess.run([os.path.join(self.name_ve, "Scripts", "activate")])
        
    def leave_ve(self):
        subprocess.run(["deactivate"])