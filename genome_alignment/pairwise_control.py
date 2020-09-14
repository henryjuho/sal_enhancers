import subprocess
import time
import argparse


def no_jobs():

    js = 'squeue -u bartonhe | wc -l'
    n = int(subprocess.Popen(js, shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[0])-1

    return n


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-sh', help='shell script to run', action='append', required=True)
    args = parser.parse_args()

    list_of_arrays = args.sh

    for a in list_of_arrays:

        # bide time while jobs run
        while no_jobs() > 50:
            time.sleep(30)

        # once break out of that loop
        sub = 'sbatch --account=tuyida {}'.format(a)
        subprocess.call(sub, shell=True)

        time.sleep(120)  # need previous jobs to start running


if __name__ == '__main__':
    main()
