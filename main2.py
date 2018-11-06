import generic as gen
import containers as cont
import sim_ops_containers as simopc
import ops_containers as opsc
import file_ops as fo
import time

def main():

    description = ""
    args = gen.parse_arguments(description, ["input_file1", "input_file2", "output_directory", "extract_sequences", "clean_run"], flags = [3,4], ints=[])
    input_file1, input_file2, output_directory, extract_sequences, clean_run = args.input_file1, args.input_file2, args.output_directory, args.extract_sequences, args.clean_run

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)


    # get the sequences
    if extract_sequences:
        cont.extract_sequences(input_file1, input_file2, output_directory, clean_run = clean_run)


if __name__ == "__main__":
    main()
