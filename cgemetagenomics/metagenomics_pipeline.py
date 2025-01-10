import os
import sys
import subprocess
import csv
import shutil

from cgemetagenomics import kma
from cgemetagenomics import version

# New PDF report imports
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS
from pathlib import Path
from datetime import datetime
import base64
import shutil
import matplotlib.pyplot as plt
import matplotlib
from io import BytesIO
import numpy as np
from collections import defaultdict

def metagenomics_pipeline(args):
    print("Starting the metagenomics pipeline...")

    # Check if output folder already exists
    output_dir = '/var/lib/cge/results/{}'.format(args.name)
    if os.path.exists(output_dir):
        sys.exit(
            f"Error: Output directory '{output_dir}' already exists. Please choose a different name or delete the existing directory.")

    if args.db_dir is None:
        if not os.path.exists('/var/lib/cge/database/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /var/lib/cge/database/cge_db')
        else:
            args.db_dir = '/var/lib/cge/database/cge_db'
            print(f"Using CGE database directory: {args.db_dir}")

    if args.output is None:
        args.output = '/var/lib/cge/results/{}'.format(args.name)

    # Create output directory
    print(f"Creating output directory: {args.output}")
    os.system('mkdir -p ' + args.output)

    # Create output directory
    # Copy the input FASTQ file to the output directory
    print(f"Copying input FASTQ file to the output directory: {args.input} -> {args.output}")
    shutil.copy(args.input, args.output)

    # Load pathogen species
    print("Loading pathogen species...")
    species = load_pathogen_species(args.db_dir + '/pathogen_strains.list')

    # Run KMA for bacteria alignment
    print("Running KMA for bacteria alignment...")
    kma.KMARunner(args.input,
                  args.output + "/bacteria_alignment",
                  args.db_dir + '/bac_db/bac_db',
                  "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8 -ef").run()

    bacterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")
    

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            e_coli_depth = find_max_depth_for_escherichia_coli(bacterial_results)
            print(f"Escherichia coli detected with depth {e_coli_depth}. Running virulence finder...")
            kma.KMARunner(args.input,
                          args.output + "/virulence",
                          args.db_dir + '/virulence_db/virulence_db',
                          "-ont -md {} -mem_mode -t 8".format(e_coli_depth / 2)).run()
            break
    
    # Run KMA for AMR
    print("Running KMA for AMR...")
    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    # Generate report
    print("Creating refined report...")
    report = create_refined_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results, species, args.name)
    with open(args.output + '/report.txt', 'w') as report_file:
        report_file.write(report)
    
    ############################  NEW PDF REPORT  ###############################
    
    file_stats = get_file_stats(args.input, args.output)  
    
    print("Creating PDF report...")
    pdf_path = create_pdf_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results, species, args.name, file_stats)
    if pdf_path:
        print(f"PDF report generated and stored in {pdf_path}")
        
    #############################################################################

    print("Metagenomics pipeline completed successfully. Report generated and stored in " + args.output + '/report.txt')
    return 'metagenomics_pipeline'


def merge_fastq_files_unix(source_directory, output_name):
    """
    Merge all fastq.gz files in the given directory using Unix commands and save the output with the specified name in the home directory.

    Args:
    source_directory (str): Path to the directory containing fastq.gz files.
    output_name (str): Name for the output file.
    """
    # Home directory path
    home_directory = os.path.expanduser('~')

    # Output file path with the specified name
    output_file = os.path.join(home_directory, f'{output_name}.fastq.gz')

    # Creating the Unix command for concatenation
    cmd = f'cat {source_directory}/*.fastq.gz > {output_file}'

    # Executing the command
    subprocess.run(cmd, shell=True, check=True)

    print(f"All files merged into {output_file}")

def create_refined_report(phenotype_file, output, bacterial_results, species, name):
    gene_data = read_tab_separated_file(phenotype_file)
    amr_results = read_tab_separated_file(output + '/amr.res')
    e_coli_found = any(extract_species(result.get('#Template', '')) == 'Escherichia coli' for result in bacterial_results)

    pathogens = []
    non_pathogens = []
    for result in bacterial_results:
        species_name = extract_species(result.get('#Template', ''))
        if species_name in species:
            pathogens.append(result)
        else:
            non_pathogens.append(result)

    phenotypes = set()
    for amr_result in amr_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result.get('#Template'):
                phenotypes.update(gene['Phenotype'].split(','))

    report = "Analysis Report: {}\n".format(name)
    report += "=" * 60 + "\n"

    # AMR Results Section
    report += "Antimicrobial Resistance (AMR) Findings:\n"
    report += "-" * 60 + "\n"
    for result in amr_results:
        report += f"Template: {result['#Template']}\n"
        report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"

    # Pathogens Found Section
    report += "Identified Pathogens:\n"
    report += "-" * 60 + "\n"
    for pathogen in pathogens:
        report += f"Template: {pathogen['#Template']}\n"
        report += f"Depth: {pathogen['Depth'].strip()}, Coverage: {pathogen['Template_Coverage'].strip()}, Identity: {pathogen['Template_Identity'].strip()}, Length: {pathogen['Template_length'].strip()}\n\n"

    # Non-Pathogens Section
    report += "Non-Pathogenic Bacteria Detected:\n"
    report += "-" * 60 + "\n"
    for non_pathogen in non_pathogens:
        report += f"Template: {non_pathogen['#Template']}\n"
        report += f"Depth: {non_pathogen['Depth'].strip()}, Coverage: {non_pathogen['Template_Coverage'].strip()}, Identity: {non_pathogen['Template_Identity'].strip()}, Length: {non_pathogen['Template_length'].strip()}\n\n"

    # Expected Phenotypes Based on AMR Genes Section
    report += "Expected Phenotypes Based on AMR Genes:\n"
    report += "-" * 60 + "\n"
    if phenotypes:
        for phenotype in sorted(phenotypes):
            report += f"• {phenotype.strip()}\n"
    else:
        report += "No phenotypes expected based on AMR genes.\n"
    report += "\n"

    # Virulence Factors for Escherichia coli Section
    if e_coli_found:
        virulence_results = read_tab_separated_file(output + '/virulence.res')
        report += "Virulence Factors for Escherichia coli:\n"
        report += "-" * 60 + "\n"
        for result in virulence_results:
            report += f"Template: {result['#Template']}\n"
            report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
    else:
        report += "No virulence factors analysis for Escherichia coli.\n"

    return report

def find_max_depth_for_escherichia_coli(bacterial_results):
    max_depth = 0.0

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            # Convert depth to a float and compare with the current max
            depth = float(hit['Depth'].strip())
            if depth > max_depth:
                max_depth = depth

    return max_depth

def load_tsv(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)


def load_pathogen_species(strain_file):
    strains = []
    with open(strain_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            species = line[1].split(' ')[0] + ' ' + line[1].split(' ')[1]
            strains.append(species)
    return strains


def read_tab_separated_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)

def extract_species(template_str):
    return ' '.join(template_str.split()[1:3])









############################  NEW PDF REPORT  ###############################


def check_report_dependencies():
    """Check if required packages for PDF generation are installed"""
    try:
        import jinja2
        import weasyprint
        return True
    except ImportError as e:
        print(f"Missing required packages for PDF report generation: {str(e)}")
        print("Please install required packages: pip install jinja2 weasyprint")
        return False

def prepare_template_data(bacterial_results, amr_results, phenotypes, pathogens, name, gene_data, file_stats):
    """Prepare comprehensive data for the template"""
    # Generate plots
    coverage_plot = generate_coverage_plot(bacterial_results)
    circle_plot = generate_amr_circle_plot(amr_results, gene_data)
    
    # Calculate AMR class statistics
    class_stats = defaultdict(lambda: {'genes': [], 'depths': [], 'coverages': []})
    
    for result in amr_results:
        gene_id = result['#Template']
        depth = float(result['Depth'].strip())
        coverage = float(result['Template_Coverage'].strip())
        
        for gene in gene_data:
            if gene['Gene_accession no.'] == gene_id:
                classes = [c.strip() for c in gene['Class'].split(',')]
                for drug_class in classes:
                    class_stats[drug_class]['genes'].append(gene_id)
                    class_stats[drug_class]['depths'].append(depth)
                    class_stats[drug_class]['coverages'].append(coverage)
    
    # Format class statistics
    amr_class_stats = []
    for class_name, stats in sorted(class_stats.items()):
        amr_class_stats.append({
            'name': class_name,
            'gene_count': len(set(stats['genes'])),
            'avg_depth': f"{sum(stats['depths']) / len(stats['depths']):.1f}",
            'avg_coverage': f"{sum(stats['coverages']) / len(stats['coverages']):.1f}"
        })
    
    return {
        'name': name,
        'date': datetime.now().strftime("%B %d, %Y"),
        'file_stats': file_stats if file_stats else {
            'file_name': 'Unknown',
            'file_size': 'N/A',
            'read_count': 'N/A',
            'kma_version': 'N/A'
        },
        'pathogen_count': len(pathogens),
        'amr_count': len(amr_results),
        'class_count': len(class_stats),  # Changed from phenotype_count
        'version': version.__version__,
        'coverage_plot': coverage_plot,
        'circle_plot': circle_plot,
        'pathogens': [
            {
                '#Template': p['#Template'],
                'Template_Coverage': p['Template_Coverage'].strip(),
                'Depth': p['Depth'].strip(),
                'Template_Identity': p['Template_Identity'].strip()
            } 
            for p in pathogens
        ],
        'amr_results': sort_amr_results(amr_results, gene_data),
        'amr_class_stats': amr_class_stats  # New field replacing grouped_phenotypes
    }



def sort_amr_results(results, gene_data):
    """Sort AMR results by Class and add quality indicators"""
    gene_to_class = {gene['Gene_accession no.']: gene['Class'] for gene in gene_data}
    
    enriched_results = []
    for result in results:
        template = result['#Template']
        identity = float(result['Template_Identity'].strip())
        coverage = float(result['Template_Coverage'].strip())

        # Exact matching with debug output
        if identity >= 100.0 and coverage >= 100.0:
            match_quality = 'perfect-match'
        elif identity >= 100.0 or coverage >= 100.0:
            match_quality = 'good-match'
        else:
            match_quality = 'partial-match'
            
        enriched_results.append({
            '#Template': template,
            'Class': gene_to_class.get(template, 'Unknown'),
            'Template_Identity': f"{identity:.1f}",
            'Template_Coverage': f"{coverage:.1f}",
            'Depth': f"{float(result['Depth'].strip()):,.1f}",
            'match_quality': match_quality
        })
    
    return sorted(enriched_results, key=lambda x: (x['Class'], -float(x['Template_Identity'])))

def create_pdf_report(phenotype_file, output, bacterial_results, species, name, file_stats):
    """Generate a PDF report alongside the text report"""
    # Check dependencies first
    if not check_report_dependencies():
        print("Skipping PDF report generation due to missing dependencies")
        return None
    
    # Get the same data as the text report
    gene_data = read_tab_separated_file(phenotype_file)
    amr_results = read_tab_separated_file(output + '/amr.res')
    
    # Setup paths
    package_dir = Path(__file__).parent
    assets_dir = package_dir / 'assets'
    logo_path = assets_dir / 'dtu_logo.png'
    
    # Verify logo exists and load it
    if not logo_path.exists():
        print(f"Warning: Logo file not found at {logo_path}")
        return None
    
    # Load and encode logo
    with open(logo_path, 'rb') as f:
        logo_data = f.read()
        logo_base64 = base64.b64encode(logo_data).decode('utf-8')
    
    # Prepare pathogen lists
    pathogens = []
    non_pathogens = []
    for result in bacterial_results:
        species_name = extract_species(result.get('#Template', ''))
        if species_name in species:
            pathogens.append(result)
        else:
            non_pathogens.append(result)
    
    # Get phenotypes
    phenotypes = set()
    for amr_result in amr_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result.get('#Template'):
                phenotypes.update(gene['Phenotype'].split(','))
    
    # Create reports directory if it doesn't exist
    reports_dir = Path(output) / 'reports'
    reports_dir.mkdir(exist_ok=True)
    
    # Prepare template data
    template_data = prepare_template_data(
        bacterial_results=bacterial_results,
        amr_results=amr_results,
        phenotypes=phenotypes,
        pathogens=pathogens,
        name=name,
        gene_data=gene_data,
        file_stats=file_stats
    )
    
    # Add logo data to template data
    template_data['logo_data_url'] = f'data:image/png;base64,{logo_base64}'
        
    # Setup Jinja2 environment
    env = Environment(
        loader=FileSystemLoader(package_dir / 'templates'),
        autoescape=True
    )
    template = env.get_template('report.html')
    
    # Render HTML
    html_content = template.render(**template_data)
    
    # Create PDF
    pdf_path = Path(output) / f"{name}_report.pdf"
    HTML(string=html_content).write_pdf(
        pdf_path,
        stylesheets=[CSS(assets_dir / 'style.css')]
    )
    
    return pdf_path

def get_extended_dtu_colors(n_colors_needed, max_colors=15):
    """Generate an extended color palette based on DTU colors"""
    import colorsys
    
    # Base DTU colors
    base_colors = [
        '#980001',  # red
        '#2F3EEA',  # blue
        '#1FD082',  # green
        '#030F4F',  # navy
        '#79238E',  # purple
        '#F6D04D',  # yellow
        '#FC7634',  # orange
        '#F7BBB1',  # light red
        '#E83F48',  # bright red
        '#008835',  # dark green
    ]
    
    def create_color_variation(hex_color, variation_index):
        """Create a variation of a color by adjusting its lightness and saturation"""
        rgb = tuple(int(hex_color.lstrip('#')[i:i+2], 16)/255 for i in (0, 2, 4))
        h, l, s = colorsys.rgb_to_hls(*rgb)
        
        if variation_index % 2 == 0:
            new_l = min(1, l * 1.3)
            new_s = max(0, s * 0.8)
        else:
            new_l = max(0, l * 0.7)
            new_s = min(1, s * 1.2)
            
        rgb = colorsys.hls_to_rgb(h, new_l, new_s)
        return '#{:02x}{:02x}{:02x}'.format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        )
    
    if n_colors_needed <= len(base_colors):
        return base_colors[:n_colors_needed]
    
    if n_colors_needed > max_colors:
        n_colors_needed = max_colors
    
    extended_colors = base_colors.copy()
    variation_index = 0
    while len(extended_colors) < n_colors_needed:
        base_color = base_colors[variation_index % len(base_colors)]
        new_color = create_color_variation(base_color, len(extended_colors))
        if new_color not in extended_colors:
            extended_colors.append(new_color)
        variation_index += 1
    
    return extended_colors

def get_plot_colors(n_categories, include_other=False, max_colors=15):
    """Get colors for plotting, handling the 'Other' category appropriately"""
    OTHERS_COLOR = '#DADADA'
    
    if n_categories > max_colors:
        n_colors_needed = max_colors - 1 if include_other else max_colors
    else:
        n_colors_needed = n_categories - 1 if include_other else n_categories
    
    colors = get_extended_dtu_colors(n_colors_needed, max_colors)
    if include_other:
        colors.append(OTHERS_COLOR)
    return colors

def generate_coverage_plot(bacterial_results):
    """Generate an enhanced coverage plot with extended DTU colors"""
    
    matplotlib.rcParams['figure.dpi'] = 300
    matplotlib.rcParams['figure.figsize'] = [4, 3]
    
    # Process top 5 bacterial results
    species_coverage = {}
    for result in bacterial_results[:5]:
        species = ' '.join(result['#Template'].split()[1:3])
        coverage = float(result['Template_Coverage'].strip())
        species_coverage[species] = coverage
    
    # Create figure
    fig = plt.figure(figsize=(4, 3), dpi=300)
    ax = fig.add_subplot(111)
    
    # Get colors for the number of species we have
    colors = get_plot_colors(len(species_coverage))
    
    # Create bar plot
    bars = ax.bar(species_coverage.keys(), species_coverage.values(),
                 color=colors)
    
    # Style improvements
    ax.set_ylabel('Coverage (%)', fontsize=8)
    ax.set_title('Top Bacterial Species Coverage', fontsize=10, pad=10)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout(pad=0.5)
    
    # Save plot
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight',
                pad_inches=0.1, facecolor='white', edgecolor='none')
    plt.close()
    
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')

def generate_amr_circle_plot(amr_results, gene_data):
    """Generate a donut chart showing AMR distribution by drug class"""
    
    MAX_CATEGORIES = 15  # Maximum categories before using "Other"
    THRESHOLD = 0.05  # 5% threshold for significance
    
    # Calculate total depth for each drug class
    class_depths = defaultdict(float)
    for result in amr_results:
        gene_id = result['#Template']
        depth = float(result['Depth'].strip())
        
        for gene in gene_data:
            if gene['Gene_accession no.'] == gene_id:
                classes = [c.strip() for c in gene['Class'].split(',')]
                for drug_class in classes:
                    class_depths[drug_class] += depth
    
    # Sort and prepare data
    sorted_data = sorted(class_depths.items(), key=lambda x: x[1], reverse=True)
    total_depth = sum(d for _, d in sorted_data)
    
    # Prepare main data and handle "Other"
    main_data = []
    main_labels = []
    other_sum = 0
    
    for drug_class, depth in sorted_data:
        percentage = depth/total_depth
        if percentage >= THRESHOLD and len(main_labels) < MAX_CATEGORIES - 1:
            main_data.append(depth)
            main_labels.append(drug_class)
        else:
            other_sum += depth
    
    if other_sum > 0:
        main_data.append(other_sum)
        main_labels.append('Other')
    
    # Get colors
    colors = get_plot_colors(len(main_labels), include_other='Other' in main_labels)
    
    # Create figure
    plt.figure(figsize=(8, 8), dpi=300, facecolor='white')
    ax = plt.gca()
    
    # Create donut chart
    wedges, texts, autotexts = ax.pie(
        main_data,
        labels=main_labels,
        colors=colors,
        autopct=lambda pct: f'{pct:.1f}%' if pct >= 5 else '',
        pctdistance=0.85,
        wedgeprops=dict(width=0.5)
    )
    
    plt.setp(autotexts, size=9, weight='bold', color='white')
    plt.setp(texts, size=10)
    
    plt.title('AMR Gene Distribution\nby Drug Class', 
              pad=20, 
              size=14, 
              weight='bold')
    
    note = "Percentages represent relative abundance based on read depth"
    plt.figtext(0.5, 0.02, note,
                ha='center',
                fontsize=8,
                style='italic')
    
    plt.tight_layout()
    
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight',
                facecolor='white')
    plt.close()
    
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def get_file_stats(input_file, output_dir):
    """Get file stats using quick commands and mapstat info.
    
    If some values cannot be found (e.g., mapstat not found or missing entries),
    only those values are 'N/A'. The absence of one file does not revert all 
    values to defaults.
    """
    file_name = 'Unknown'
    file_size = 'N/A'
    fragment_count = None
    kma_version = None
    
    # Attempt to get file_name and file_size from input file
    if os.path.exists(input_file):
        file_name = os.path.basename(input_file)
        try:
            file_size_cmd = f"ls -lh {input_file} | cut -d' ' -f5"
            file_size = subprocess.check_output(file_size_cmd, shell=True).decode().strip()
            if not file_size:
                file_size = 'N/A'
        except:
            file_size = 'N/A'
    
    # Attempt to get fragment_count and kma_version from mapstat
    mapstat_file = os.path.join(output_dir, "bacteria_alignment.mapstat")
    if os.path.exists(mapstat_file):
        try:
            with open(mapstat_file, 'r') as f:
                for line in f:
                    if line.startswith("## fragmentCount"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            try:
                                fragment_count = int(parts[1])
                            except ValueError:
                                fragment_count = None
                    elif line.startswith("## version"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            kma_version = parts[1]
        except:
            # If parsing fails for some reason, just continue with what we have
            pass

    return {
        'file_name': file_name if file_name else 'Unknown',
        'file_size': file_size if file_size else 'N/A',
        'read_count': f"{fragment_count:,}" if fragment_count is not None else 'N/A',
        'kma_version': kma_version if kma_version else 'N/A'
    }