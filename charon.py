"""Ferry the souls of pickled output files across the river SSHtyx off
of tara and back to the local machine."""

import os, subprocess, shlex, string

i = 0
while(True):
    files = os.listdir("upstream")
    pickles = [f for f in files if f.endswith("pickle")]
    if pickles:
        archive_name = "archive%s.tar" % i
        tar_command = "tar -czf %s %s" % (achive_name, " ".join(pickles))
        subprocess.call(shlex.split(tar_command))
        scp_command = "scp "
        i += 1
