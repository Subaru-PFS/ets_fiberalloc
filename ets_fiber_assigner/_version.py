import subprocess

__version__ = "3.2.0"

try:
    tmp = subprocess.run(["git", "describe", "--dirty"], capture_output=True, text=True)
    __gitversion__ = tmp.stdout.rstrip()
except:
    __gitversion__ = "unknown"

def recursive_version_info():
    import ics.cobraOps
    import pfs.utils
    print(f"ets_fiber_assigner {__version__} ({__gitversion__}), depending on {{")
    try:
        ics.cobraOps.recursive_version_info()
    except:
        try:
            print(f"ics.cobraOps {ics.cobraOps.__version__}")
        except:
            print("ics.cobraOps (unknown)")
    try:
        pfs.utils.recursive_version_info()
    except:
        try:
            print(f"pfs.utils {pfs.utils.__version__}")
        except:
            print("pfs.utils (unknown)")
    print("}")
