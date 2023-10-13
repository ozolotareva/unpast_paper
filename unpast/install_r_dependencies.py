import os
import sys
import subprocess

def check_r_installation():
    try:
        subprocess.run(["R", "--version"], check=True, capture_output=True)
        return True
    except FileNotFoundError:
        print("R is not installed on this system. Skipping R library installation.")
        return False

def install_bioc_manager():
    try:
        subprocess.run(["R", "-e", "install.packages('BiocManager')"], check=True)
    except subprocess.CalledProcessError:
        print("An error occurred while installing 'BiocManager'.")
        sys.exit(1)

def install_r_library(library):
    try:
        subprocess.run(["R", "-e", f"BiocManager::install('{library}')"], check=True)
        print(f"R library '{library}' has been installed successfully.")
    except subprocess.CalledProcessError:
        print(f"An error occurred while installing the R library '{library}'.")
        sys.exit(1)

if __name__ == "__main__":
    if check_r_installation():
        install_bioc_manager()
        for library in ['WGCNA', 'limma']:
            install_r_library(library)

