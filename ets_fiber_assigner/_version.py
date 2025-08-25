import subprocess

__version__ = "3.0"

try:
    tmp = subprocess.run(["git", "describe", "--dirty"], capture_output=True, text=True)
    gitversion = tmp.stdout
except:
    gitversion = "unknown.gitversion"

__version__ = __version__ + "+" + gitversion
