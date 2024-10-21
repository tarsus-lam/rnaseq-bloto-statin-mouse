import os
import numpy as np
import matplotlib.pyplot as plt

# Parse alignment statistics from STAR log file
def parse_star_log(log_file):
    alignment_stats = {}
    with open(log_file, 'r') as f:
        for line in f:
            if '|' in line:
                parts = line.strip().split('|')
                alignment_stats[parts[0].strip()] = parts[1].strip().replace(',', '')
    return alignment_stats

# Calculate summary statistics
def calculate_summary(alignment_stats_list):
    total_mapped_reads = sum([int(stats['Uniquely mapped reads number']) for stats in alignment_stats_list])
    total_multi_mapped_reads = sum([int(stats['Number of reads mapped to multiple loci']) for stats in alignment_stats_list])
    total_reads = sum([int(stats['Number of input reads']) for stats in alignment_stats_list])

    percent_mapped = (total_mapped_reads / total_reads) * 100
    percent_unmapped = ((total_reads - total_mapped_reads - total_multi_mapped_reads) / total_reads) * 100
    percent_multi_mapped = (total_multi_mapped_reads / total_reads) * 100
    
    return percent_mapped, percent_unmapped, percent_multi_mapped

# Plot alignment summary
def plot_alignment_summary(percent_mapped, percent_unmapped, percent_multi_mapped, log_directory):
    labels = ['Mapped', 'Unmapped', 'Multi-mapped']
    sizes = [percent_mapped, percent_unmapped, percent_multi_mapped]
    colors = ['lightblue', 'lightcoral', 'lightskyblue']
    explode = (0.1, 0, 0)  # explode 1st slice

    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title('Average Alignment Summary')
    
    # Check if "AllignmentStats" directory exists, if not, create it
    if not os.path.exists('AlignmentStats'):
        os.makedirs('AlignmentStats')
    
    # Save the plot in "AlignmentStats" directory
    plt.savefig(os.path.join(outdir, 'alignment_summary.png'))
    plt.show()

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python script.py log_directory output_directory")
        sys.exit(1)

    log_directory = sys.argv[1]
    outdir = sys.argv[2]

    # List all files in the directory
    log_files = [f for f in os.listdir(log_directory) if f.endswith("_Log.final.out")]

    # Parse alignment statistics from all log files
    alignment_stats_list = [parse_star_log(os.path.join(log_directory, file)) for file in log_files]

    # Calculate average summary statistics
    percent_mapped, percent_unmapped, percent_multi_mapped = calculate_summary(alignment_stats_list)

    # Plot average alignment summary and save plot
    plot_alignment_summary(percent_mapped, percent_unmapped, percent_multi_mapped, outdir)

