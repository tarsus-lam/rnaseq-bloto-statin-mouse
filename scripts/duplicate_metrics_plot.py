import os
import sys
import matplotlib.pyplot as plt

def parse_metrics_file(metrics_file):
    """
    Parse the metrics file and extract relevant metrics.
    """
    metrics = {}
    with open(metrics_file, 'r') as file:
        for line in file:
            if line.startswith('## METRICS CLASS'):
                metrics_class = line.split()[-1]
                metrics[metrics_class] = {}
                # Parse metric name from the first line after '## METRICS CLASS'
                metric_name_line = next(file).strip().split('\t')
                # Parse metric value from the second line
                metric_value_line = next(file).strip().split('\t')
                # Store metric name and value in the dictionary
                for name, value in zip(metric_name_line, metric_value_line):
                    metrics[metrics_class][name] = value
                break  # Stop parsing after finding metrics class
    return metrics

def plot_summary(metrics_directory, output_directory):
    """
    Plot summary metrics for all files in the directory.
    """
    metrics_files = [f for f in os.listdir(metrics_directory) if f.endswith('_marked_metrics.txt')]
    total_reads = []
    duplicate_reads = []
    sample_names = []
    for file in metrics_files:
        metrics = parse_metrics_file(os.path.join(metrics_directory, file))
        total = float(metrics['picard.sam.DuplicationMetrics']['READ_PAIRS_EXAMINED']) + float(metrics['picard.sam.DuplicationMetrics']['UNPAIRED_READS_EXAMINED'])
        dup = float(metrics['picard.sam.DuplicationMetrics']['READ_PAIR_DUPLICATES']) + float(metrics['picard.sam.DuplicationMetrics']['UNPAIRED_READ_DUPLICATES'])
        total_reads.append(total)
        duplicate_reads.append(dup)
        sample_names.append(os.path.splitext(file)[0])
    
    # Calculate the percentage of duplicate reads for each sample
    percentage_duplicates = [(dup / total) * 100 for total, dup in zip(total_reads, duplicate_reads)]
    
    # Sort samples based on their percentage of duplicates
    sorted_indices = sorted(range(len(percentage_duplicates)), key=lambda k: percentage_duplicates[k])
    sample_names = [sample_names[i] for i in sorted_indices]
    total_reads = [total_reads[i] for i in sorted_indices]
    duplicate_reads = [duplicate_reads[i] for i in sorted_indices]
    percentage_duplicates = [percentage_duplicates[i] for i in sorted_indices]

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax1.bar(range(len(sample_names)), total_reads, color='tab:blue', label='Total Reads')
    ax1.bar(range(len(sample_names)), duplicate_reads, color='tab:orange', label='Duplicate Reads')
    ax2.plot(range(len(sample_names)), percentage_duplicates, color='tab:green', marker='o', label='% Duplicate Reads')

    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Number of Reads', color='black')
    ax2.set_ylabel('% Duplicate Reads', color='black')
    ax1.set_title('Duplicate Metrics Summary')

    ax1.set_xticks(range(len(sample_names)))
    ax1.set_xticklabels(sample_names, rotation=45, ha='right')

    # Add legend for both bar plots and the line plot
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper left')

    fig.tight_layout()
    output_path = os.path.join(output_directory, 'duplicate_metrics_summary_plot.png')
    plt.savefig(output_path)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python duplicate_metrics_plot.py <bam_directory> <output_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]

    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Plot summary
    plot_summary(input_directory, output_directory)
