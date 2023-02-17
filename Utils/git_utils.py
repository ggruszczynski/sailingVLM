from subprocess import CalledProcessError
import subprocess
import os

def get_git_revision_hash() -> str:
    try:
        output = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip())
    except subprocess.CalledProcessError:
        return None
    except FileNotFoundError:
        print("No git repository found.", "ERROR")
        # log("No git repository found.", "ERROR")
        return None

    return output

def get_git_revision_short_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()

def get_git_branch(path=None):
    if path is None:
        path = os.path.curdir
    command = 'git rev-parse --abbrev-ref HEAD'.split()
    try:
        branch = subprocess.Popen(command, stdout=subprocess.PIPE, cwd=path).stdout.read()
        output = str(branch.strip().decode('utf-8'))
    except subprocess.CalledProcessError:
        return None
    except FileNotFoundError:
        print("No git repository found.", "ERROR")
        # log("No git repository found.", "ERROR")
        return None
    return output
