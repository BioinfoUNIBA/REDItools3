import sys


def check_dead(processes):
    """
    Look through processes to determine if any have died unexpectedly.

    If any process has an exit code of 1, this method will terminate all other
    processes and then exit with code 1.

    Parameters:
        processes (list): Processes to check
    """
    for proc in processes:
        if proc.exitcode == 1:
            for to_kill in processes:
                to_kill.kill()
            sys.stderr.write('[ERROR] Killing job\n')
            sys.exit(1)
