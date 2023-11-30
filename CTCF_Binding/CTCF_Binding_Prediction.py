import pyft
import tqdm

def read_bed_file(bed_file):
    """
    Read regions from a BED file.

    :param bed_file: Path to the BED file.
    :return: Dictionary with regions (key: region name, value: (chromosome, start, end, strand)).
    """
    regions = {}
    with open(bed_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # Skip invalid lines

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[5] if len(parts) > 5 else '+'  # Default to '+' strand if not specified
            region_name = f"{chrom}:{start}-{end}"
            regions[region_name] = (chrom, start, end, strand)
    return regions

def analyze_ctcf_binding_nucleosomes(bam_file, regions):
    """
    Analyze the binding status of CTCF Zinc finger regions in fiberseq data considering the DNA strand.
    
    :param bam_file: Path to the BAM file.
    :param regions: Dictionary with regions of interest (key: region name, value: (chromosome, start, end, strand)).
    :return: Dictionary with binding status and nucleosome count for each region.
    """
    fiberbam = pyft.Fiberbam(bam_file)
    binding_status = {}

    for region_name, (chrom, start, end, strand) in regions.items():
        bound = 0
        unbound = 0
        nucleosome_count = 0

        for fiber in tqdm.tqdm(fiberbam.fetch(chrom, start, end)):
            # Skip fibers that don't match the strand
            if fiber.strand != strand:
                continue

            # Check for nucleosomes
            if is_nucleosome(fiber):
                nucleosome_count += 1
                continue

            # Check for m6a methylation marks in the region
            if has_m6a_methylation(fiber):
                unbound += 1
            else:
                bound += 1

        binding_status[region_name] = {'bound': bound, 'unbound': unbound, 'nucleosomes': nucleosome_count}

    return binding_status

def has_m6a_methylation(fiber):
    """
    Check if a fiber has m6a methylation marks.
    
    :param fiber: Fiberdata object.
    :return: True if m6a methylation marks are present, False otherwise.
    """
    for start in fiber.m6a.reference_starts:
        if start >= fiber.start and start <= fiber.end:
            return True
    return False

def is_nucleosome(fiber):
    """
    Determine if a fiber represents a nucleosome.
    
    :param fiber: Fiberdata object.
    :return: True if the fiber is identified as a nucleosome, False otherwise.
    """
    # Implement nucleosome identification logic here
    return any(start >= fiber.start and start <= fiber.end for start in fiber.nuc.reference_starts)

# Path to the BED file and BAM file
bed_file = 'GM12878_CTCF_overlap.bed'
bam_file = 'PS0015i.norm.bam'

# Read regions from the BED file
regions = read_bed_file(bed_file)

# Analyze binding and nucleosomes
binding_analysis = analyze_ctcf_binding_nucleosomes(bam_file, regions)

# Print the results
print(binding_analysis)
