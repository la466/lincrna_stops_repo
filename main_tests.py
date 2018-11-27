import generic as gen
import containers as cont
import sim_ops_containers as simopc
import ops_containers as opsc
import file_ops as fo
import time

def main():

    arguments = ["input_file1", "input_file2", "output_directory", "stop_density"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [3], ints=[])
    input_file1, input_file2, output_directory, stop_density = args.input_file1, args.input_file2, args.output_directory, args.stop_density

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    # get the stop density
    if stop_density:
        # input_file1 = fasta containing sequences, input_file2 = bed file containing families
        sequo.get_stop_density(input_file1, families_file = input_file2)


if __name__ == "__main__":
    main()
