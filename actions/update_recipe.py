
import sys
import os
import subprocess

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"kcombu_bss source is in {srcdir}")

condadir = os.path.join(srcdir, "recipes", "kcombu_bss")

print(f"conda recipe in {condadir}")

# Store the name of the recipe and template YAML files.
recipe = os.path.join(condadir, "meta.yaml")
template = os.path.join(condadir, "template.yaml")

if not os.path.exists(template):
    print(f"No template recipe in {template}?")
    sys.exit(-1)


def run_cmd(cmd):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()

gitdir = os.path.join(srcdir, ".git")

print(f"git dir is {gitdir}")

# Get the Sire version. (Latest tag.)
kcombu_version = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} describe --tags --abbrev=0")
print(kcombu_version)

# Get the build number. (Number of commits since last tag.)
kcombu_build = len(run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} log --oneline {kcombu_version}..").split("\n"))
print(kcombu_build)

# Get the branch.
kcombu_branch = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} rev-parse --abbrev-ref HEAD")
print(kcombu_branch)

lines = open(template, "r").readlines()

with open(recipe, "w") as FILE:
    for line in lines:
        line = line.replace("KCOMBU_VERSION", kcombu_version)
        line = line.replace("KCOMBU_BUILD", str(kcombu_build))
        line = line.replace("KCOMBU_BRANCH", kcombu_branch)

        FILE.write(line)
