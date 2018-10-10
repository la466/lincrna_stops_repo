import generic as gen
import sim_ops as simo
import seq_ops as seqo
import os

def sim_orf_length(seq_fasta, required_simulations, output_file, parallel=True, seeds=None, seq_seeds=None):
    """
    Simulation to look at the lengths of ORFs in sequences
    compared with simulated null sequences

    Args:
        seq_fasta (str): path to fasta file containing sequences
        required_simulations (int): number of simulations to do
        output_file (str): path to output file
        parallel (bool): opt whether to use multiprocessing
        seeds (list): seed numbers for the simulations
        seq_seeds (list): a list of lists containing seeds for testing
    """

    names, seqs_list = gen.read_fasta(seq_fasta)

    # get a dictionary of unique sequences
    unique_seqs = []
    seqs = {}
    for i, name in enumerate(names):
        if seqs_list[i] not in unique_seqs:
            unique_seqs.append(seqs_list[i])
            seqs[name] = seqs_list[i]

    # get the dicnucleotide and nucleotide content of the sequences
    dinucleotide_content = seqo.get_dinucleotide_content(unique_seqs)
    nucleotide_content = seqo.get_nucleotide_content(unique_seqs)

    # get the longest orfs and gc content for the real sequences
    longest_orfs = seqo.get_longest_orfs(seqs)
    gc_contents = {}
    for seq in seqs:
        gc_contents[seq] = seqo.calc_seq_gc(seqs[seq])

    # create a temporary output directory
    temp_dir = "temp_data"
    gen.create_output_directories(temp_dir)

    # run the simulations
    simulations = list(range(required_simulations))
    sim_args = [seqs, dinucleotide_content, nucleotide_content, temp_dir, seeds, seq_seeds]

    # run the simulations
    if parallel:
        # add foo to argument list for parallelisation
        sim_args.insert(0, "foo")
        workers = os.cpu_count() - 2
        processes = gen.run_in_parallel(simulations, sim_args, simo.sim_orf_lengths, workers = workers)
        outputs = []
        for process in processes:
            outputs.extend(process.get())
    else:
        outputs = simo.sim_orf_lengths(simulations, *sim_args)

    # get temp filelist so we can order simulants for tests
    temp_filelist = {}
    for file in outputs:
        simulation_no = file.split('.')[-2]
        temp_filelist[simulation_no] = file

    # write to file
    with open(output_file, "w") as outfile:
        # write the header
        outfile.write("sim_no,{0}\n".format(",".join(name for name in sorted(seqs))))
        # write the gc content
        outfile.write("gc,{0}\n".format(",".join(gen.stringify([gc_contents[id] for id in sorted(gc_contents)]))))
        # write the real line
        outfile.write("real,{0}\n".format(",".join(gen.stringify([longest_orfs[id] for id in sorted(longest_orfs)]))))
        # write the simulants
        for sim in sorted(temp_filelist):
            # get the simulation number
            data = seqo.fasta_to_dict(temp_filelist[sim])
            outfile.write("{0},{1}\n".format("sim_{0}".format(sim), ",".join([data[id] for id in sorted(data)])))

    # remove the temp files
    gen.remove_files(outputs)
