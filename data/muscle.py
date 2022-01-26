import os
import subprocess

path_to_muscle_executable = '/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6'
script_dir = os.path.dirname(os.path.realpath(__file__))  # path to this file
data_dir = script_dir + '/spike_proteins/'  # relative path of datasets


def muscle_align_fasta_file(infilename,
                            outfilename,
                            data_dir = data_dir,
                            path_to_muscle_executable ='/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle'):

    bash_command = f'{path_to_muscle_executable} -in {data_dir + infilename} -out {data_dir + outfilename} -log {data_dir + "muscle.log"}'
    print(bash_command)
    subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)


if __name__ == "__main__":
    muscle_align_fasta_file(infilename='1_in_500_cleaned.fasta', outfilename='1_in_500_cleaned_aligned.afa')
