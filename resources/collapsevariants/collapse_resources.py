import dxpy
import subprocess


# This function runs a command on an instance, either with or without calling the docker instance we downloaded
# By default, commands are not run via Docker, but can be changed by setting is_docker = True
# Also, by default, standard out is not saved, but can be modified with the 'stdout_file' parameter.
# print_cmd is for internal debugging purposes when testing new code
def run_cmd(cmd: str, is_docker: bool = False, stdout_file: str = None, print_cmd = False) -> None:

    # -v here mounts a local directory on an instance (in this case the home dir) to a directory internal to the
    # Docker instance named /test/. This allows us to run commands on files stored on the AWS instance within Docker.
    # This looks slightly different from other versions of this command I have written as I needed to write a custom
    # R script to run STAAR. That means we have multiple mounts here to enable this code to find the script.
    if is_docker:
        cmd = "docker run " \
              "-v /home/dnanexus:/test " \
              "-v /usr/bin/:/prog " \
              "egardner413/mrcepid-burdentesting " + cmd

    if print_cmd:
        print(cmd)

    # Standard python calling external commands protocol
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout_file is not None:
        with open(stdout_file, 'w') as stdout_writer:
            stdout_writer.write(stdout.decode('utf-8'))
        stdout_writer.close()

    # If the command doesn't work, print the error stream and close the AWS instance out with 'dxpy.AppError'
    if proc.returncode != 0:
        print("The following cmd failed:")
        print(cmd)
        print("STDOUT follows\n")
        print(stdout.decode('utf-8'))
        print("STDERR follows\n")
        print(stderr.decode('utf-8'))
        raise dxpy.AppError("Failed to run properly...")


# This is to generate a global CHROMOSOMES variable for parallelisation
def get_chromosomes(is_snp_tar: bool = False, is_gene_tar: bool = False):

    if is_snp_tar:
        chromosomes = list(['SNP'])
    elif is_gene_tar:
        chromosomes = list(['GENE'])
    else:
        chromosomes = list(range(1,23)) # Is 0 based on the right coordinate...? (So does 1..22)
        chromosomes.extend(['X'])
        chromosomes = list(map(str, chromosomes))

    return chromosomes
