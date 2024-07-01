"""
Rerun Systre with verbose output

Calls the Systre command, writing its full, raw stdout and stderr instead of
parsing out the RCSR topology or other error states.  This code can sometimes
be a useful diagnostic after running bin/sbu.exe or Python/run_mofid.py
"""

import sys
from mofid.id_constructor import SYSTRE_CMD_LIST, DEFAULT_SYSTRE_CGD
if sys.version_info[0] < 3:
    try:
        import subprocess32 as subprocess
    except:
        raise AssertionError('You must install subprocess32 if running Python2')
else:
    import subprocess

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) > 1:
        raise SyntaxError('Usage: python rerun_systre.py optional_path.cgd')
    if len(args) == 1:
        cgd_path = args[0]
    else:
        cgd_path = DEFAULT_SYSTRE_CGD
    
    # Run Systre without a timeout
    cmd_list = SYSTRE_CMD_LIST + [cgd_path]
    java_run = subprocess.run(cmd_list, universal_newlines=True)
    # Writing to stdout and stderr by default since they're not captured with
    # a pipe here, unlike in run_mofid.py:runcmd
    
    # Alternatively, manually re-forward the output to stdout and stderr
    #print(java_run.stdout)
    #sys.stderr.write(java_run.stderr)
