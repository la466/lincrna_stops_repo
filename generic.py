'''
Author: Rosina Savisaar and Liam Abrahams.
Module that contains generic utility functions that make life a bit easier.
'''

import argparse
import csv
import ftplib
import itertools as it
import multiprocessing
import numpy as np
import os
import random
import re
import shutil
import subprocess
import time

def create_dir(path, strict=None):
    '''
    Create new directory
    strict: if true, delete directory and make new
    '''
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    else:
        shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)

def create_directory(path, strict=None):
    '''
    Create set of directories for a given path
    strict: if true, delete directory and make new
    '''
    path_splits = path.split('/')
    new_path = []
    for i, split in enumerate(path_splits):
        new_path.append(split)
        create_dir("/".join(new_path), strict=strict)

def copy_file(src, dest):
    '''
    Copy a file from one diretory to another
    '''
    shutil.copyfile(src, dest)

def flatten(structured_list):
    '''
    Flatten a structured list.
    '''
    flat_list = list(it.chain(*structured_list))
    return(flat_list)

def get_time(start_time):
    '''
    Print out how many minutes have passed since start_time.
    '''
    current = time.time()
    spent = round((current - start_time)/60, 2)
    print("{0} minutes.\n".format(spent))

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def list_to_dict(input_list, index1, index2, as_list = False, uniquify = False, floatify = False):
    '''
    Convert the input_list into a dictionary, with the index1th element of each sublist as the key and the index2th element as the value.
    '''
    if as_list and floatify:
        print("_as_list_ and _floatify_ can't both be True!")
        raise Exception
    output_dict = {}
    for i in input_list:
        if not as_list:
            if floatify:
                output_dict[i[index1]] = float(i[index2])
            else:
                output_dict[i[index1]] = i[index2]
        else:
            if i[index1] not in output_dict:
                output_dict[i[index1]] = []
            output_dict[i[index1]].append(i[index2])
    if as_list and uniquify:
        output_dict = {i: sorted(list(set(output_dict[i]))) for i in output_dict}
    return(output_dict)

def motif_to_regex(motifs):
    '''
    Convert a string into a lookahead regex where only the first base
    is matched and the rest is in the lookahead.
    '''
    regex = [re.compile("".join([i[0],"(?=",i[1:],")"])) for i in motifs]
    return(regex)

def parse_arguments(description, arguments, floats = None, flags = None, ints = None):
    '''
    Use argparse to parse a set of input arguments from the command line.
    '''
    if not floats:
        floats = []
    if not flags:
        flags = []
    if not ints:
        ints = []
    parser = argparse.ArgumentParser(description = description)
    for pos, argument in enumerate(arguments):
        if pos in flags:
            parser.add_argument("--{0}".format(argument), action = "store_true", help = argument)
        else:
            if pos in floats:
                curr_type = float
            elif pos in ints:
                curr_type = int
            else:
                curr_type = str
            parser.add_argument(argument, type = curr_type, help = argument)
    args = parser.parse_args()
    return(args)

def read_fasta(input_file):
    '''
    Given a fasta file return a first lists containing the sequence identifiers and a second list containing teh sequences (in the same order).
    '''
    file_to_read = open(input_file)
    input_lines = file_to_read.readlines()
    file_to_read.close()
    input_lines = [i.rstrip("\n") for i in input_lines]
    names = [i.lstrip(">") for i in input_lines if i[0] == ">"]
    sequences = [i for i in input_lines if i[0] != ">"]
    if len(sequences) != len(names):
        print("Problem extracting data from fasta file!")
        print(len(sequences))
        print(len(names))
        raise Exception
    if len(sequences) == 0:
        print("No sequences were extracted!")
        raise Exception
    return(names, sequences)

def read_many_fields(input_file, delimiter):
    '''
    Read a csv/tsv/... into a list of lists with each sublist corresponding to one line.
    '''
    file_to_read = open(input_file)
    try:
        field_reader = csv.reader(file_to_read, delimiter = delimiter)
        lines = []
        for i in field_reader:
            lines.append(i)
        file_to_read.close()
        return(lines)
    except:
        print("Problem reading file...")
        return [["Problem reading file"]]

def remove_directory(dir):
    '''
    Remove directory
    '''
    if os.path.exists(dir):
        shutil.rmtree(dir)

def remove_file(file_name):
    '''
    Remove a file, if it exists.
    '''
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass

def reverse_complement(base):
    '''
    Reverse complement a base.
    '''
    reverse_comps = {"A": "T","C": "G","G": "C","T": "A"}
    return(reverse_comps[base])


def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example, a list of gene names)
    args: a list of arguments to the function. Put in "foo" in place of the argument you are parallelizing over.
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for elem in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = elem
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)


def run_process(arguments, return_string = True, input_to_pipe = None, return_error = False, file_for_input = None, file_for_output = None, univ_nl = True, shell = False):
    '''
    Run a command on the command line. Supply command as a list of strings.
    EX: run_process(["cat", "hello!"], file_for_output = "hello.txt")
    '''
    if file_for_input:
        input_file = open(file_for_input)
        stdin_src = input_file
    else:
        stdin_src = subprocess.PIPE
    if file_for_output:
        output_file = open(file_for_output, "w")
        stdout_dest = output_file
    else:
        stdout_dest = subprocess.PIPE
    arguments = [str(i) for i in arguments]
    if shell:
        arguments = " ".join(arguments)
    process = subprocess.Popen(arguments, shell = shell, stdout = stdout_dest, stderr = subprocess.PIPE,
                               stdin = stdin_src, universal_newlines = univ_nl)
    if input_to_pipe:
        stdout, stderr = process.communicate(input_to_pipe)
    else:
        stdout, stderr = process.communicate()
    if file_for_input:
        input_file.close()
    if file_for_output:
        output_file.close()
    return_code = process.poll()
    if return_code != 0:
        print("Process failed!")
        print(" ".join(arguments))
        print(stderr)
        return("error")
    #if the process returns bytes but you want to get a string back.
    if return_string and type(stdout) == bytes:
        stdout = stdout.decode("utf-8")
    if return_error:
        return(stderr)
    else:
        return(stdout)


def write_to_fasta(names, seq, fasta_name):
    '''
    Write a set of sequence identifiers and sequences to fasta file.
    '''
    with open(fasta_name, "w") as file:
        for i in range(len(names)):
            file.write(">{0}\n".format(names[i]))
            file.write("{0}\n".format(seq[i]))

def stringify(item):
    if isinstance(item, list):
        return [str(i) for i in item]
    else:
        return str(item)
