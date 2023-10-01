
import os

class CreateCommand:
    def __init__(self, module_dir, path_to_mixcr, paired_end_sequencing, experiment, files, java_heap_size, threads = None):
        self.java_heap_size = java_heap_size
        self.module_dir = module_dir
        self.path_to_mixcr = path_to_mixcr
        self.threads = threads
        self.paired_end_sequencing = paired_end_sequencing
        if paired_end_sequencing:
            self.files = [os.path.normpath(j) for i in files for j in i]
            basename = os.path.basename(os.path.splitext(files[0][0])[0])
        else:
            self.files = files
            basename = os.path.basename(os.path.splitext(files[0])[0])
        self.basename = basename
        self.result = os.path.join(self.module_dir,
                                   "temp",
                                   self.basename + ".vdjca")
        self.alignment_path = os.path.normpath(os.path.join(self.module_dir,
                                                            "my_experiments",
                                                            experiment,
                                                            "alignment_reports",
                                                            self.basename + "_AlignmentReport.txt"))
        self.assembly_path = os.path.normpath(os.path.join(self.module_dir,
                                                           "my_experiments",
                                                           experiment,
                                                           "assembly_reports",
                                                           self.basename + "_AssemblyReport.txt"))
        self.clones = os.path.join(self.module_dir,
                                   "my_experiments",
                                   experiment,
                                   "clones_result",
                                   self.basename + "_clones.clns")
        self.clones_base = os.path.dirname(self.clones)
        self.table_tsv = os.path.join(module_dir,
                                     "my_experiments",
                                     experiment,
                                     "tables_mixcr",
                                     self.basename + ".tsv")
        self.mixcr_plots_path = os.path.normpath(os.path.join(self.module_dir,
                                                              "my_experiments",
                                                              experiment,
                                                              "mixcr_plots"))

    def read_aligned_reads(self):
        with open(self.alignment_path, "r") as f:
            alignment = f.readlines()
        for i in alignment:
            if "Successfully aligned reads" in i:
                quality = i.replace("Successfully aligned reads: ", "")
        align_quality = quality.split(" ")[0]
        percent_quality = quality.split(" ")[1]
        align_quality = int(align_quality)
        return align_quality, percent_quality

    def read_clones(self):
        with open(self.assembly_path, "r") as f:
            alignment = f.readlines()
        for i in alignment:
            if "Reads used in clonotypes, percent of total" in i:
                quality = i.replace("Reads used in clonotypes, percent of total: ", "")
        align_quality = quality.split(" ")[0]
        percent_quality = quality.split(" ")[1]
        assembly_quality = int(align_quality)
        return assembly_quality, percent_quality

    def get_free_memory_linux(self):
        with open("/proc/meminfo", "r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("MemAvailable"):
                return int(line.split()[1]) // 1024  # Convert from KB to MB

    def get_free_memory_windows(self):
        try:
            import psutil
            return psutil.virtual_memory().available // (1024 ** 2)  # Convert from Bytes to MB
        except ImportError:
            print("Please install psutil library: pip install psutil")


    def get_free_memory_mac(self):
        cmd = "vm_stat | grep 'Pages free' | awk '{print $3}'"
        pages_free = int(os.popen(cmd).read().replace(".", ""))
        # Typically, a page in macOS is 4096 bytes
        return (pages_free * 4096) // (1024 ** 2)  # Convert from Bytes to MB



    def create_parser(self):
     #   try:
      #      free_memory = self.get_free_memory()
      #      java_heap_size = int(free_memory / 4)  # Using half of the available RAM as an example
       #     print(f"25 % of currently available memory: {java_heap_size}\n")
        #except:
         #   print("automatic detection of memory failed. Processing continues with 1000MB RAM.")
          #  java_heap_size = 1000
        commands = []
        commands.extend(["java", f"-Xms{self.java_heap_size}M", "-jar"]) # enable change of para
        commands.extend([self.path_to_mixcr])
        return commands


    def prepare_align(self, method):
        align_commands = self.create_parser()
        align_commands.extend(["align"])
        print(method)
        align_commands.extend(["--preset", method])
        align_commands.extend(self.files)
        align_commands.extend([self.result])
        align_commands.extend(["--no-warnings"])
        align_commands.extend(["--force-overwrite"])
        if self.threads != None:
            align_commands.extend(["--threads", self.threads])
        align_commands.extend(["--report", self.alignment_path])
        return align_commands

    def prepare_assembly(self):

        assembly_commands = self.create_parser()
        assembly_commands.extend(["assemble"])
       # assembly_commands.extend(["-OseparateByC=true", "-OseparateByV=true","-OseparateByJ=true"])
        assembly_commands.extend([self.result])
        assembly_commands.extend([self.clones])
        assembly_commands.extend(["--force-overwrite"])
        if self.threads != None:
            assembly_commands.extend(["--threads", self.threads])
        assembly_commands.extend(["--report", self.assembly_path])

        return assembly_commands


    def prepare_clones(self):
        clones_commands = self.create_parser()
        clones_commands.extend(["exportClones"])
        clones_commands.extend(["-c IGH"])
        clones_commands.extend([self.clones])
        clones_commands.extend(["--force-overwrite"])
        if self.threads != None:
            clones_commands.extend(["--threads", self.threads])
        clones_commands.extend([self.table_tsv])
        return clones_commands

