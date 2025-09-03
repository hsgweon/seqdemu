#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import multiprocessing
import os, shutil
import tempfile
import progressbar
from multiprocessing import Manager


parser = argparse.ArgumentParser("seqdemu: A No-Nonsense Nanopore Demultiplexer.")
parser.add_argument("-i",
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED] Input gzipped FASTQ file.", 
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED] Output file prefix.",
                    required = True)
parser.add_argument("-b",
                    action = "store",
                    dest = "barcodes",
                    metavar = "barcodes",
                    help = "[REQUIRED] Barcode file in TSV format (fwd_primer, rev_primer, sample_name).",
                    required = True)
parser.add_argument("-m",
                    action = "store",
                    dest = "mismatch",
                    metavar = "mismatch",
                    help = "[REQUIRED] Number of exact mismatches allowed for barcodes.",
                    required = True)
parser.add_argument("-t",
                    action = "store",
                    dest = "threads",
                    metavar = "threads",
                    help = "[REQUIRED] Number of threads to use.",
                    required = True)
parser.add_argument("-f", "--format",
                    action="store",
                    dest="format",
                    default="both",
                    choices=['fasta', 'fastq', 'both'],
                    help="Output format: 'fasta', 'fastq', or 'both'. Default: 'both'")
options = parser.parse_args()


def boyer_moore_search_with_exact_mismatches(text, pattern, exact_mismatches):
    n = len(text)
    m = len(pattern)
    results = []

    # Preprocess pattern for Boyer-Moore
    bad_character = {}
    for i in range(m - 1):
        bad_character[ord(pattern[i])] = m - i - 1

    # Search with Boyer-Moore
    i = 0
    while i <= n - m:
        j = m - 1
        mismatches = 0
        while j >= 0:
            if pattern[j] != text[i + j]:
                mismatches += 1
                if mismatches > exact_mismatches:
                    break
            j -= 1

        if j == -1 and mismatches == exact_mismatches:
            results.append((i, i + m - 1))
            i += 1
        else:
            i += max(bad_character.get(ord(text[i + m - 1]), m), m - j - 1)

    return results

def process_chunk(args):
    chunk, progress_counter, lock, temp_file_prefix_cpu = args

    # Determine which formats we need to write based on user option
    formats_to_write = []
    if options.format in ['fasta', 'both']:
        formats_to_write.append('fasta')
    if options.format in ['fastq', 'both']:
        formats_to_write.append('fastq')

    # Create a dictionary to hold file handles for each format and trim type
    temp_files = {fmt: {} for fmt in formats_to_write}
    for fmt in formats_to_write:
        temp_files[fmt]['full'] = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_full_{i}.{fmt}"), "w") for i in range(len(list_barcodes))]
        temp_files[fmt]['barc'] = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_barc_{i}.{fmt}"), "w") for i in range(len(list_barcodes))]
        temp_files[fmt]['noba'] = [open(os.path.join(temp_dir, f"{temp_file_prefix_cpu}_noba_{i}.{fmt}"), "w") for i in range(len(list_barcodes))]

    exact_mismatches = int(options.mismatch)

    for records in chunk:
        for record in records:
            sequence = str(record.seq).upper()

            for i in range(len(list_barcodes)):
                # Strand 0
                strand = 0
                barcode_F_0 = list_barcodes[i][strand][0]
                barcode_RC_0 = list_barcodes[i][strand][1]
                result_F_0 = boyer_moore_search_with_exact_mismatches(sequence, barcode_F_0, exact_mismatches)
                result_RC_0 = boyer_moore_search_with_exact_mismatches(sequence, barcode_RC_0, exact_mismatches)

                # Strand 1
                strand = 1
                barcode_F_1 = list_barcodes[i][strand][0]
                barcode_RC_1 = list_barcodes[i][strand][1]
                result_F_1 = boyer_moore_search_with_exact_mismatches(sequence, barcode_F_1, exact_mismatches)
                result_RC_1 = boyer_moore_search_with_exact_mismatches(sequence, barcode_RC_1, exact_mismatches)
                
                # Check if there is just one forward and one reverse primer
                STRAND_0_VALID = (len(result_F_0) == 1 and len(result_RC_0) == 1) and (result_F_0[0][0] < result_RC_0[0][0])
                STRAND_1_VALID = (len(result_F_1) == 1 and len(result_RC_1) == 1) and (result_F_1[0][0] < result_RC_1[0][0])
                
                # Define records to write if a valid match is found
                record_to_write = None
                if STRAND_0_VALID and not STRAND_1_VALID:
                    start_barc, end_barc = result_F_0[0][0], result_RC_0[0][1] + 1
                    start_noba, end_noba = result_F_0[0][1] + 1, result_RC_0[0][0]
                    record_to_write = True
                elif not STRAND_0_VALID and STRAND_1_VALID:
                    start_barc, end_barc = result_F_1[0][0], result_RC_1[0][1] + 1
                    start_noba, end_noba = result_F_1[0][1] + 1, result_RC_1[0][0]
                    record_to_write = True

                # If a valid hit was found, create sliced records and write them
                if record_to_write:
                    # Slicing a SeqRecord object slices sequence and quality scores together
                    barc_record = record[start_barc:end_barc]
                    noba_record = record[start_noba:end_noba]
                    
                    # Set the description to match the original for consistency
                    barc_record.description = record.description
                    noba_record.description = record.description

                    # Write the records in the requested formats using SeqIO.write
                    for fmt in formats_to_write:
                        # For FASTQ, specify the quality format
                        format_name = 'fastq-sanger' if fmt == 'fastq' else 'fasta'
                        
                        SeqIO.write(record, temp_files[fmt]['full'][i], format_name)
                        SeqIO.write(barc_record, temp_files[fmt]['barc'][i], format_name)
                        SeqIO.write(noba_record, temp_files[fmt]['noba'][i], format_name)

            with lock:
                progress_counter.value += 1
                pbar.update(progress_counter.value)

    # Close all temporary file handles
    for fmt in formats_to_write:
        for trim_type in ['full', 'barc', 'noba']:
            for f_handle in temp_files[fmt][trim_type]:
                f_handle.close()


def divide_file_chunks(infile, chunk_size):
    chunks = []
    # Input is always FASTQ, so this part is unchanged
    with gzip.open(infile, "rt") as f:
        records = list(SeqIO.parse(f, "fastq"))
        for i in range(0, len(records), chunk_size):
            chunks.append(records[i:i + chunk_size])
    return chunks


if __name__ == "__main__":
    
    list_barcodes = []
    list_samplenames = []

    # Count number of sequences in FASTQ
    with gzip.open(options.infile, 'rb') as f:
        # A bit of a rough way to count, but it works for the progress bar
        num_lines = sum(1 for line in f)
    number_of_sequences = int(num_lines / 4)

    manager = Manager()
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Barcodes - deals with Forward-RC(Reverse) AND Reverse-RC(Forward)
    with open(options.barcodes, "r") as barcodes:
        for line in barcodes:
            parts = line.rstrip().split("\t")
            list_barcodes.append([(parts[0], str(Seq(parts[1]).reverse_complement())), (parts[1], str(Seq(parts[0]).reverse_complement()))])
            list_samplenames.append(parts[2])

    print("Calculating how to divide FASTQ file into chunks and assigning to multiple CPUs...")

    num_cpus = int(options.threads)
    chunk_size = 1000 # Number of sequences in a chunk
    chunks = divide_file_chunks(options.infile, chunk_size)

    temp_dir = "temp_seqdemu"
    shutil.rmtree(temp_dir, ignore_errors=True)
    os.makedirs(temp_dir)

    pbar = progressbar.ProgressBar(max_value=number_of_sequences).start()

    pool_args = [(chunks[i::num_cpus], progress_counter, lock, f"cpu{i}") for i in range(num_cpus)]

    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.map(process_chunk, pool_args)

    pbar.finish()
    print("Demultiplexing complete. Merging temporary files...")

    
    # Determine which formats we need to process based on user option
    formats_to_process = []
    if options.format in ['fasta', 'both']:
        formats_to_process.append('fasta')
    if options.format in ['fastq', 'both']:
        formats_to_process.append('fastq')

    # Loop through each format and merge the corresponding files
    for fmt in formats_to_process:
        print(f"Generating .{fmt} output files...")
        outfiles = {}
        # Open final output files with the correct extension
        for name in list_samplenames:
            outfiles[name] = {
                "full": open(f"{options.outfile}_{name}_full.{fmt}", "w"),
                "barc": open(f"{options.outfile}_{name}_barc.{fmt}", "w"),
                "noba": open(f"{options.outfile}_{name}_noba.{fmt}", "w")
            }

        # Concatenate temporary files into final output files
        for trim in ["full", "barc", "noba"]:
            for i in range(len(list_samplenames)):
                for cpu in range(num_cpus):
                    temp_filepath = os.path.join(temp_dir, f"cpu{cpu}_{trim}_{i}.{fmt}")
                    if os.path.exists(temp_filepath):
                        with open(temp_filepath, "r") as temp_file:
                            # Efficiently copy file contents
                            shutil.copyfileobj(temp_file, outfiles[list_samplenames[i]][trim])
        
        # Close all final output files for this format
        for name in list_samplenames:
            for trim in ["full", "barc", "noba"]:
                outfiles[name][trim].close()
                

    shutil.rmtree(temp_dir, ignore_errors=True)
    
    print("All done.")
