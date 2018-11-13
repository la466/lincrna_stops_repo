import generic as gen
import containers as cont
import sim_ops_containers as simopc
import ops_containers as opsc
import file_ops as fo
import time

def main():

    arguments = ["input_file1", "input_file2", "input_file3", "input_file4", "input_file5", "output_directory", "extract_sequences", "clean_run"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [6,7], ints=[])
    input_file1, input_file2, input_file3, input_file4, input_file5, output_directory, extract_sequences, clean_run = args.input_file1, args.input_file2, args.input_file3, args.input_file4, args.input_file5, args.output_directory, args.extract_sequences, args.clean_run

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    # get the sequences
    if extract_sequences:
        # input_file1 = gtf genome 1, input_file2 = genome fasta 1, input_file3 = gtf genome 2, input_file4 = genome fasta 2, input_file5 = orthlogs file
        # cont.extract_clean_sequences(input_file1, input_file2, input_file3, input_file4, input_file5, output_directory, clean_run = clean_run)
        cont.filter_families(input_file1)

if __name__ == "__main__":
    main()
