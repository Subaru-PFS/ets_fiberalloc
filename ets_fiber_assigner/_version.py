import subprocess

__version__ = "3.1.0"

try:
    tmp = subprocess.run(["git", "describe", "--dirty"], capture_output=True, text=True)
    __gitversion__ = tmp.stdout.rstrip()
except:
    __gitversion__ = "unknown"
