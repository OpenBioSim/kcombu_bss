
import os
import sys
import glob

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"Source is in {srcdir}\n")

# Get the anaconda token to authorise uploads
if "ANACONDA_TOKEN" in os.environ:
    conda_token = os.environ["ANACONDA_TOKEN"]
else:
    conda_token = "TEST"

# get the build directory
if "BUILD_DIR" in os.environ:
    conda_bld = os.environ["BUILD_DIR"]
else:
    conda_bld = os.path.join("..", "build")

print(f"conda_bld = {conda_bld}")

# Find the packages to upload
kcombu_pkg = glob.glob(os.path.join(conda_bld, "*-*", "kcombu*.tar.bz2"))

if len(kcombu_pkg) == 0:
    print("No sire packages to upload?")
    sys.exit(-1)

packages = kcombu_pkg

print(f"Uploading packages:")
print(" * ", "\n *  ".join(packages))

packages = " ".join(packages)

def run_cmd(cmd):
    import subprocess
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()

print(f"\nLabelling with 'main' and 'dev'.")
label = "--label main --label dev"

# Upload the packages to the michellab channel on Anaconda Cloud.
cmd = f"anaconda --token {conda_token} upload --user michellab {label} --force {packages}"

print(f"\nUpload command:\n\n{cmd}\n")

if conda_token == "TEST":
    print("Not uploading as the ANACONDA_TOKEN is not set!")
    sys.exit(-1)

output = run_cmd(cmd)

print(output)

print("Package uploaded!")
