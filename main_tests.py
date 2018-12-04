import generic as gen
import containers as cont
import sim_ops_containers as simopc
import main_tests_ops as mto
import ops_containers as opsc
import file_ops as fo
import time

def main():

    arguments = ["input_directory", "output_directory", "compare_stop_density", "compare_stop_density_limit_frames", "coding_exons", "generate_gc_controls_exons", "generate_gc_controls_introns", "generate_dint_exon_controls", "cds_ds", "cds_density_nd"]

    description = ""
    args = gen.parse_arguments(description, arguments, flags = [2,3,4,5,6,7,8,9], ints=[])
    input_directory, output_directory, compare_stop_density, compare_stop_density_limit_frames, coding_exons, generate_gc_controls_exons, generate_gc_controls_introns, generate_dint_exon_controls, cds_ds, cds_density_nd = args.input_directory, args.output_directory, args.compare_stop_density, args.compare_stop_density_limit_frames, args.coding_exons, args.generate_gc_controls_exons, args.generate_gc_controls_introns, args.generate_dint_exon_controls, args.cds_ds, args.cds_density_nd

    # set a start time
    start = time.time()

    # create the output_directory if it doenst already exist
    gen.create_output_directories(output_directory)

    exons_fasta = "{0}/genome_sequences/human/human.cds.clean_coding_exons.fasta".format(input_directory)
    cds_fasta = "{0}/genome_sequences/human/human.cds.clean.fasta".format(input_directory)
    introns_fasta = "{0}/genome_sequences/human/human.clean_introns.fasta".format(input_directory)
    non_transcribed_fasta = "{0}/genome_sequences/human/human.non_transcribed.fasta".format(input_directory)
    families_file = "{0}/genome_sequences/human/human.cds.families.bed".format(input_directory)


    gc_control_exon_output_directory = "{0}/clean_exon_gc_controls".format(output_directory)
    if generate_gc_controls_exons:
        simopc.generate_gc_controls(exons_fasta, non_transcribed_fasta, gc_control_exon_output_directory)

    gc_control_intron_output_directory = "{0}/clean_intron_gc_controls".format(output_directory)
    if generate_gc_controls_introns:
        simopc.generate_gc_controls(introns_fasta, non_transcribed_fasta, gc_control_intron_output_directory)

    dint_control_cds_output_directory = "{0}/clean_cds_dint_controls".format(output_directory)
    if generate_dint_exon_controls:
        simopc.generate_dint_controls(cds_fasta, dint_control_cds_output_directory)

    # get the stop density
    output_file = "{0}/compare_exon_intron_density.csv".format(output_directory)
    if compare_stop_density:
        mto.compare_stop_density(exons_fasta, introns_fasta, output_file, families_file = families_file)

    


    # if coding_exons:
    #     mto.coding_exons(input_file1, input_file2, output_directory)

    # cds_fasta = "{0}/genome_sequences/human/human.cds.clean.fasta".format(input_directory)
    # ortholog_cds_fasta = "{0}/genome_sequences/macaque/macaque.cds.quality_filtered.step1.fasta".format(input_directory)
    # ortholog_transcript_links = "{0}/genome_sequences/human.macaque.conservation_filtering.step3.bed".format(input_directory)
    # if cds_ds:
    #     mto.position_ds(cds_fasta, ortholog_cds_fasta, ortholog_transcript_links, output_directory)
    #
    #
    # gc_controls_zip = "{0}/clean_coding_exons_gc_controls.zip".format(input_directory)
    #
    # if cds_density_nd:
    #     mto.cds_density_nd(exons_fasta, families_file, gc_controls_zip, output_directory)

if __name__ == "__main__":
    main()
