import sys
import os
import pysam
import matplotlib.pyplot as plt
import concurrent.futures

# Function to calculate read coverage from a BAM file
def calculate_coverage(bam_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    coverage = [0] * samfile.lengths[0]  # Initialize coverage list with zeros
    for pileupcolumn in samfile.pileup():
        coverage[pileupcolumn.pos] = pileupcolumn.n  # Update coverage at each position
    samfile.close()
    return coverage

# Function to merge coverage data from multiple BAM files
def merge_coverage(bam_files, num_threads):
    max_len = max(pysam.AlignmentFile(bam_file, "rb").lengths[0] for bam_file in bam_files)
    merged_coverage = [0] * max_len  # Initialize merged coverage list with zeros
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit coverage calculation tasks to the executor
        coverage_futures = [executor.submit(calculate_coverage, bam_file) for bam_file in bam_files]
        # Merge coverage from all BAM files
        for future in concurrent.futures.as_completed(coverage_futures):
            coverage = future.result()
            for i in range(len(coverage)):
                merged_coverage[i] += coverage[i]  # Sum coverage from all BAM files
    return merged_coverage

# Function to plot aggregated read coverage
def plot_coverage(coverage, output_dir):
    plt.plot(range(1, len(coverage) + 1), coverage)
    plt.xlabel('Genomic Position')
    plt.ylabel('Coverage')
    plt.title('Aggregated Read Coverage Plot')
    output_path = os.path.join(output_dir, 'read_coverage_plot.png')
    plt.savefig(output_path)
    print(f"Plot saved to: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <bam_directory> <output_directory> <num_threads>")
        sys.exit(1)

    # Input directory containing BAM files
    bam_directory = sys.argv[1]

    # Desired output directory to save the plot
    output_directory = sys.argv[2]

    # Number of threads for parallel processing
    num_threads = int(sys.argv[3])

    # List all BAM files in the input directory
    bam_files = ['STARAlignmentsPass2/N2_S34_L004_Aligned.sortedByCoord.out.bam', 'STARAlignmentsPass2/N15_S47_L004_Aligned.sortedByCoord.out.bam']
    #bam_files = [os.path.join(bam_directory, f) for f in os.listdir(bam_directory) if f.endswith('.bam')]

    # Merge coverage data from all BAM files
    merged_coverage = merge_coverage(bam_files, num_threads)

    # Plot aggregated read coverage and save the plot
    plot_coverage(merged_coverage, output_directory)

