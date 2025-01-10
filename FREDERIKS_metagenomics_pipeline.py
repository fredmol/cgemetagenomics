import os
import sys
import subprocess
import csv
import shutil
import os
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.piecharts import Pie
from reportlab.graphics.charts.legends import Legend
from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Frame, PageTemplate, Flowable
from reportlab.lib.colors import HexColor, Color
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.graphics.shapes import Circle, Drawing, String
from reportlab.lib.units import mm
from reportlab.platypus import PageBreak, Frame, PageTemplate, Flowable
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime


from cgemetagenomics import kma

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
                  "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8").run()

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

    # Generate reports
    print("Creating refined report...")
    report = create_refined_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results, species, args.name)
    with open(args.output + '/report.txt', 'w') as report_file:
        report_file.write(report)
         
    print("Creating PDF report...")
    pdf_path = create_beautiful_pdf_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results, species, args.name)

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


def create_beautiful_pdf_report(phenotype_file, output, bacterial_results, species, name):
    """Creates a beautiful PDF report with visualizations"""
    # Get the same data as create_refined_report
    gene_data = read_tab_separated_file(phenotype_file)
    amr_results = read_tab_separated_file(output + '/amr.res')
    e_coli_found = any(extract_species(result.get('#Template', '')) == 'Escherichia coli' 
                      for result in bacterial_results)
    
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
                
    virulence_results = None
    if e_coli_found:
        virulence_results = read_tab_separated_file(output + '/virulence.res')

    # Generate PDF
    report_generator = BeautifulReportGenerator(output, name)
    return report_generator.create_beautiful_report(
        bacterial_results=bacterial_results,
        amr_results=amr_results,
        phenotypes=phenotypes,
        pathogens=pathogens,
        non_pathogens=non_pathogens,
        virulence_results=virulence_results
    )


class BeautifulReportGenerator:
    def __init__(self, output_dir, name):
        self.output_dir = output_dir
        self.name = name
        self.pdf_path = os.path.join(output_dir, f"{name}_beautiful_report.pdf")
        self.styles = getSampleStyleSheet()
        self.colors = {
            'primary': HexColor('#1B4F72'),
            'secondary': HexColor('#2E86C1'),
            'accent': HexColor('#3498DB'),
            'text': HexColor('#2C3E50'),
            'light': HexColor('#ECF0F1')
        }
        self.create_custom_styles()

    def create_custom_styles(self):
        """Create magazine-style custom styles"""
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            fontName='Helvetica-Bold',
            fontSize=32,
            spaceAfter=30,
            textColor=self.colors['primary'],
            leading=36
        ))
        
        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            fontName='Helvetica-Bold',
            fontSize=20,
            spaceAfter=15,
            textColor=self.colors['secondary'],
            leading=24,
            borderPadding=(10, 0, 10, 0),
            borderWidth=0,
            borderColor=self.colors['light'],
            backColor=self.colors['light']
        ))
        
        self.styles.add(ParagraphStyle(
            name='CustomBodyText',
            fontName='Helvetica',
            fontSize=11,
            leading=16,
            textColor=self.colors['text']
        ))

    class ResultCircle(Flowable):
        def __init__(self, size, value, label, color):
            Flowable.__init__(self)
            self.size = size
            self.value = value
            self.label = label
            self.color = color
            
        def draw(self):
            # Create circle
            radius = self.size/2 - 2
            self.canv.setFillColor(self.color)
            self.canv.circle(self.size/2, self.size/2, radius, fill=1)
            
            # Add value text
            self.canv.setFillColor(colors.white)
            self.canv.setFont("Helvetica", 20)
            self.canv.drawCentredString(self.size/2, self.size/2-10, str(self.value))
            
            # Add label text
            self.canv.setFillColor(self.color)
            self.canv.setFont("Helvetica", 8)
            self.canv.drawCentredString(self.size/2, 10, self.label)

    def create_beautiful_report(self, bacterial_results, amr_results, phenotypes, 
                              pathogens, non_pathogens, virulence_results=None):
        """Generate a magazine-style PDF report"""
        doc = SimpleDocTemplate(
            self.pdf_path,
            pagesize=A4,
            rightMargin=25*mm,
            leftMargin=25*mm,
            topMargin=20*mm,
            bottomMargin=20*mm
        )
        
        story = []
        
        # Title
        story.append(Paragraph(f"Metagenomics Analysis Report", self.styles['CustomTitle']))
        story.append(Spacer(1, 10*mm))
        
        # Key metrics in table format
        data = [
            ["Key Findings"],
            [f"Pathogens Detected: {len(pathogens)}"],
            [f"AMR Genes Found: {len(amr_results)}"],
            [f"Resistance Types: {len(phenotypes)}"]
        ]
        
        table_style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), self.colors['secondary']),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 14),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('TOPPADDING', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 1, self.colors['light']),
        ])
        
        summary_table = Table(data, colWidths=[400])
        summary_table.setStyle(table_style)
        story.append(summary_table)
        story.append(Spacer(1, 10*mm))
        
        # Rest of content
        story.extend(self._create_main_content(bacterial_results, amr_results, 
                                             phenotypes, pathogens, virulence_results))
        
        doc.build(story)
        return self.pdf_path

    def _create_main_content(self, bacterial_results, amr_results, phenotypes, 
                           pathogens, virulence_results):
        """Create main report content with enhanced styling"""
        content = []
        
        # Sample info with enhanced styling
        content.append(Paragraph("Sample Overview", self.styles['SectionHeader']))
        content.append(Paragraph(self.create_pathogen_summary(pathogens), 
                               self.styles['CustomBodyText']))
        content.append(Paragraph(self.create_amr_summary(amr_results, phenotypes), 
                               self.styles['CustomBodyText']))
        content.append(Spacer(1, 10*mm))
        
        # Coverage plot with enhanced styling
        content.append(Paragraph("Bacterial Composition", self.styles['SectionHeader']))
        coverage_plot = self.generate_coverage_plot(bacterial_results)
        content.append(Image(coverage_plot, width=400, height=200))
        content.append(Spacer(1, 10*mm))
        
        # AMR profile with enhanced styling
        if phenotypes:
            content.append(Paragraph("Antimicrobial Resistance Profile", 
                                   self.styles['SectionHeader']))
            for phenotype in sorted(phenotypes):
                content.append(Paragraph(f"• {phenotype.strip()}", 
                                       self.styles['CustomBodyText']))
            content.append(Spacer(1, 10*mm))
        
        return content

    def generate_coverage_plot(self, bacterial_results):
        """Generate an enhanced coverage plot"""
        fig, ax = plt.subplots(figsize=(10, 5))
        
        species_coverage = {}
        for result in bacterial_results[:5]:
            species = ' '.join(result['#Template'].split()[1:3])
            coverage = float(result['Template_Coverage'].strip())
            species_coverage[species] = coverage
        
        # Convert HexColor to RGB tuple
        color = self.colors['secondary']
        rgb_color = (color.red, color.green, color.blue)
        
        bars = ax.bar(species_coverage.keys(), species_coverage.values(),
                     color=rgb_color)
        
        ax.set_ylabel('Coverage (%)', fontsize=12)
        ax.set_title('Top Bacterial Species Coverage', fontsize=14, pad=20)
        plt.xticks(rotation=45, ha='right')
        
        ax.grid(True, axis='y', alpha=0.3)
        plt.tight_layout()
        
        plot_path = os.path.join(self.output_dir, 'coverage_plot.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return plot_path
    
    def create_pathogen_summary(self, pathogens):
        """Create a simple summary of pathogens found"""
        if not pathogens:
            return "No pathogens were detected in this sample."
            
        pathogen_counts = len(pathogens)
        high_coverage_pathogens = sum(1 for p in pathogens 
                                    if float(p['Template_Coverage'].strip()) > 80)
        
        summary = (f"Found {pathogen_counts} potential pathogen(s). "
                  f"{high_coverage_pathogens} showed high coverage (>80%).")
        return summary

    def create_amr_summary(self, amr_results, phenotypes):
        """Create a human-readable AMR summary"""
        if not amr_results:
            return "No antimicrobial resistance genes were detected."
            
        resistance_count = len(phenotypes)
        summary = (f"Found genes indicating potential resistance to {resistance_count} "
                  f"different types of antimicrobials.")
        return summary
    
def parse_report_data(bacterial_results, amr_results, phenotypes, pathogens):
    return {
        "quick_summary": {
            "Pathogens Found": len(pathogens),
            "AMR Genes Found": len(amr_results),
            "Resistance Phenotypes": len(phenotypes)
        },
        "bacterial_results": bacterial_results,
        "amr_results": amr_results,
        "phenotypes": phenotypes,
        "pathogens": pathogens
    }

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