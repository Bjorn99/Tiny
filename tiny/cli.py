import typer
import json
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeRemainingColumn
from rich.layout import Layout
from rich.syntax import Syntax
from rich.tree import Tree
from rich import box
from rich.style import Style
from rich.rule import Rule
from typing import List, Optional
from pathlib import Path
from enum import Enum

from tiny.core.sequence import DNASequence
from tiny.core.analysis import analyze_sequences, compare_sequences
from tiny.core.advanced import AdvancedAnalysis
from tiny.core.utils import load_fasta, save_results
from tiny.core.formats import FormatHandler

app = typer.Typer(help="🧬 Tiny - DNA Sequence Analysis Tool")
console = Console()

class AlignmentMode(str, Enum):
    global_align = "global"
    local = "local"
    semi_global = "semi-global"

def create_analysis_table(title: str) -> Table:
    """Create a styled table for analysis results."""
    table = Table(
        title=title,
        box=box.ROUNDED,
        header_style="bold cyan",
        title_style="bold magenta",
        border_style="blue",
        show_edge=True
    )
    return table

def show_progress(description: str):
    """Create a progress context manager."""
    return Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(style="cyan"),
        TimeRemainingColumn(),
        console=console
    )

@app.command()
def analyze(
    sequences: List[str] = typer.Argument(None),
    input_file: Optional[Path] = typer.Option(None, "--input", "-i", 
        help="Input file (FASTA, FASTQ, GenBank, EMBL)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", 
        help="Output file for results"),
    format_info: bool = typer.Option(False, "--format-info", 
        help="Show additional format-specific information"),
    feature_limit: int = typer.Option(5, "--feature-limit", "-fl",
        help="Number of features to display per type (0 for all)"),
    feature_type: str = typer.Option(None, "--feature-type", "-ft",
        help="Filter features by type (e.g., CDS, gene, tRNA)"),
    save_features: bool = typer.Option(False, "--save-features", "-sf",
        help="Save complete feature information to a separate file")
):
    """✨ Analyze DNA sequences and display comprehensive results."""
    try:
        with show_progress("Analyzing sequences") as progress:
            task = progress.add_task("Processing...", total=100)
            
            if input_file:
                records = FormatHandler.read_file(input_file)
                progress.update(task, advance=30)
                
                table = create_analysis_table(f"Sequence Analysis Results: {input_file.name}")
                table.add_column("ID", style="cyan")
                table.add_column("Length", justify="right", style="magenta")
                table.add_column("GC Content", justify="right", style="green")
                table.add_column("Molecular Weight", justify="right", style="yellow")
                
                if format_info:
                    table.add_column("Additional Info", style="blue")
                
                progress.update(task, advance=20)
                
                results = []
                for record in records:
                    dna_seq = DNASequence(record.sequence)
                    result = analyze_sequences([dna_seq])[0]
                    results.append(result)
                    
                    row = [
                        record.id,
                        str(result.length),
                        f"{result.gc_content:.2f}%",
                        f"{result.molecular_weight:.2f}"
                    ]
                    
                    if format_info:
                        info = []
                        if record.features:
                            info.append(f"Features: {len(record.features)}")
                        if record.quality_scores:
                            avg_quality = sum(record.quality_scores) / len(record.quality_scores)
                            info.append(f"Avg Quality: {avg_quality:.2f}")
                        if record.metadata:
                            for key, value in record.metadata.items():
                                info.append(f"{key}: {value}")
                        row.append("\n".join(info))
                    
                    table.add_row(*row)
                
                progress.update(task, advance=20)
                
                console.print()
                console.rule("[bold blue]Analysis Complete", style="blue")
                console.print()
                console.print(table)
                
                if format_info and any(r.features for r in records):
                    # Group features by type
                    feature_types = {}
                    for record in records:
                        if record.features:
                            for feature in record.features:
                                feat_type = feature['type']
                                if feature_type and feat_type != feature_type:
                                    continue
                                if feat_type not in feature_types:
                                    feature_types[feat_type] = []
                                feature_types[feat_type].append((record.id, feature))
                    
                    # Show summary table
                    console.print()
                    summary_table = Table(title="Feature Type Summary", box=box.ROUNDED)
                    summary_table.add_column("Feature Type", style="cyan")
                    summary_table.add_column("Count", justify="right", style="magenta")
                    
                    for feat_type, features in sorted(feature_types.items()):
                        summary_table.add_row(feat_type, str(len(features)))
                    
                    console.print(summary_table)
                    
                    # Add detailed feature tree
                    console.print()
                    feature_tree = Tree("🧬 [bold cyan]Detailed Features Overview")
                    
                    for feat_type, features in sorted(feature_types.items()):
                        type_node = feature_tree.add(f"[blue]{feat_type} ({len(features)})")
                        # Show features based on limit
                        display_features = features if feature_limit == 0 else features[:feature_limit]
                        for record_id, feature in display_features:
                            feature_info = f"[green]{record_id} at {feature['location']}"
                            if feature['qualifiers']:
                                # Show all qualifiers
                                feature_node = type_node.add(feature_info)
                                for key, value in feature['qualifiers'].items():
                                    feature_node.add(f"[yellow]{key}: {value}")
                            else:
                                type_node.add(feature_info)
                        
                        if feature_limit > 0 and len(features) > feature_limit:
                            type_node.add(f"[yellow]... and {len(features) - feature_limit} more")
                    
                    console.print(feature_tree)

                    # Save complete feature information if requested
                    if save_features:
                        features_output = output.parent / f"{output.stem}_features.json" if output else Path("features.json")
                        with open(features_output, 'w') as f:
                            feature_data = {
                                feat_type: [
                                    {
                                        'record_id': record_id,
                                        'location': str(feature['location']),
                                        'qualifiers': feature['qualifiers']
                                    }
                                    for record_id, feature in features
                                ]
                                for feat_type, features in feature_types.items()
                            }
                            json.dump(feature_data, f, indent=2)
                        console.print(f"\n[green]Complete feature information saved to:[/green] {features_output}")

                progress.update(task, advance=30)
                
            else:
                results = analyze_sequences([DNASequence(seq) for seq in sequences])
                progress.update(task, advance=50)
                
                table = create_analysis_table("Direct Sequence Analysis")
                table.add_column("Sequence", style="cyan")
                table.add_column("Length", justify="right", style="magenta")
                table.add_column("GC Content", justify="right", style="green")
                table.add_column("Molecular Weight", justify="right", style="yellow")

                for result in results:
                    table.add_row(
                        f"{str(result.sequence)[:30]}{'...' if len(str(result.sequence)) > 30 else ''}",
                        str(result.length),
                        f"{result.gc_content:.2f}%",
                        f"{result.molecular_weight:.2f}"
                    )
                
                progress.update(task, advance=30)
                
                console.print()
                console.rule("[bold blue]Analysis Complete", style="blue")
                console.print()
                console.print(table)
            
            if output:
                save_results(results, output)
                console.print(f"\n[green]Results saved to:[/green] {output}")

    except Exception as e:
        console.print(f"\n[red bold]Error:[/red bold] {str(e)}")
        raise typer.Exit(1)

@app.command()
def compare(
    sequence1: str = typer.Argument(..., help="First DNA sequence"),
    sequence2: str = typer.Argument(..., help="Second DNA sequence"),
    input_file: Optional[Path] = typer.Option(None, "--input", "-i", 
        help="Input file containing sequences to compare")
):
    """🔍 Compare two DNA sequences for mutations and similarities."""
    try:
        with show_progress("Comparing sequences") as progress:
            task = progress.add_task("Processing...", total=100)
            
            if input_file:
                records = FormatHandler.read_file(input_file)
                if len(records) < 2:
                    raise ValueError("Input file must contain at least two sequences")
                sequence1 = records[0].sequence
                sequence2 = records[1].sequence
            
            progress.update(task, advance=50)
            comparison = compare_sequences(DNASequence(sequence1), DNASequence(sequence2))
            progress.update(task, advance=50)
            
            console.print()
            console.rule("[bold blue]Comparison Results", style="blue")
            console.print()
            
            # Create comparison results panel
            results_table = Table(show_header=False, box=box.SIMPLE)
            results_table.add_column(style="cyan")
            results_table.add_column(style="green")
            results_table.add_row("Identity:", f"{comparison.identity:.2f}%")
            results_table.add_row("Alignment Score:", f"{comparison.alignment_score}")
            
            console.print(Panel(results_table, title="[bold cyan]Sequence Comparison", border_style="blue"))
            
            if comparison.mutations:
                mutations_table = Table(title="Mutations Found", box=box.ROUNDED)
                mutations_table.add_column("Position", style="cyan", justify="right")
                mutations_table.add_column("Original", style="red")
                mutations_table.add_column("Mutated", style="green")
                
                for mut in comparison.mutations:
                    mutations_table.add_row(
                        str(mut.position),
                        mut.original,
                        mut.mutated
                    )
                
                console.print()
                console.print(mutations_table)
    
    except Exception as e:
        console.print(f"\n[red bold]Error:[/red bold] {str(e)}")
        raise typer.Exit(1)

@app.command()
def align(
    sequence1: str = typer.Argument(..., help="First DNA sequence"),
    sequence2: str = typer.Argument(..., help="Second DNA sequence"),
    mode: AlignmentMode = typer.Option(
        AlignmentMode.global_align,
        help="Alignment mode"
    ),
):
    """Perform sequence alignment with detailed visualization."""
    try:
        with show_progress("Performing alignment") as progress:
            task = progress.add_task("Processing...", total=100)
            
            progress.update(task, advance=50)
            result = AdvancedAnalysis.align_sequences(sequence1, sequence2, mode)
            progress.update(task, advance=50)
            
            console.print()
            console.rule("[bold blue]Alignment Results", style="blue")
            console.print()
            
            # Create alignment statistics panel
            stats_table = Table(show_header=False, box=box.SIMPLE)
            stats_table.add_column(style="cyan")
            stats_table.add_column(style="green")
            stats_table.add_row("Mode:", f"[blue]{mode}[/blue]")
            stats_table.add_row("Score:", f"{result.alignment_score:.2f}")
            stats_table.add_row("Identity:", f"{result.identity:.2f}%")
            stats_table.add_row("Gaps:", str(result.gaps))
            
            console.print(Panel(stats_table, title="[bold cyan]Alignment Statistics", border_style="blue"))
            
            # Show alignment visualization
            console.print()
            alignment_panel = Panel(
                f"[cyan]Sequence 1:[/cyan] {result.aligned_seq1}\n"
                f"[cyan]Sequence 2:[/cyan] {result.aligned_seq2}",
                title="[bold cyan]Alignment Visualization",
                border_style="blue"
            )
            console.print(alignment_panel)
    
    except Exception as e:
        console.print(f"\n[red bold]Error:[/red bold] {str(e)}")
        raise typer.Exit(1)

@app.command()
def find_motifs(
    sequences: List[str] = typer.Argument(None),
    input_file: Optional[Path] = typer.Option(None, "--input", "-i"),
    motif_length: int = typer.Option(6, "--length", "-l"),
    min_frequency: int = typer.Option(2, "--min-freq", "-m"),
):
    """Find common motifs in DNA sequences."""
    try:
        with show_progress("Finding motifs") as progress:
            task = progress.add_task("Processing...", total=100)
            
            if input_file:
                sequences = [str(record.sequence) for record in FormatHandler.read_file(input_file)]
            elif not sequences:
                raise ValueError("Please provide sequences or an input file")
            
            progress.update(task, advance=50)
            motifs = AdvancedAnalysis.find_motifs(sequences, motif_length, min_frequency)
            progress.update(task, advance=50)
            
            console.print()
            console.rule("[bold blue]Motif Analysis Results", style="blue")
            console.print()
            
            table = create_analysis_table("Found Motifs")
            table.add_column("Motif", style="cyan")
            table.add_column("Frequency", justify="right", style="magenta")
            table.add_column("Consensus Score", justify="right", style="green")
            table.add_column("Positions", style="yellow")

            for motif in motifs:
                table.add_row(
                    motif.motif,
                    str(motif.frequency),
                    f"{motif.consensus_score:.3f}",
                    ", ".join(map(str, motif.positions))
                )

            console.print(table)
    
    except Exception as e:
        console.print(f"\n[red bold]Error:[/red bold] {str(e)}")
        raise typer.Exit(1)

@app.command()
def find_regulatory(
    sequence: str = typer.Argument(..., help="DNA sequence to analyze"),
):
    """🧬 Find potential regulatory elements in a DNA sequence."""
    try:
        with show_progress("Analyzing regulatory elements") as progress:
            task = progress.add_task("Processing...", total=100)
            
            progress.update(task, advance=50)
            elements = AdvancedAnalysis.find_regulatory_elements(sequence)
            progress.update(task, advance=50)
            
            console.print()
            console.rule("[bold blue]Regulatory Elements Analysis", style="blue")
            console.print()
            
            found_elements = False
            for element_type, findings in elements.items():
                if findings:
                    found_elements = True
                    table = create_analysis_table(element_type.replace('_', ' ').title())
                    table.add_column("Position", style="cyan", justify="right")
                    table.add_column("Sequence", style="green")
                    
                    for pos, seq in findings:
                        table.add_row(str(pos), seq)
                    
                    console.print(table)
                    console.print()
            
            if not found_elements:
                console.print("[yellow]No regulatory elements found in the sequence.[/yellow]")
    
    except Exception as e:
        console.print(f"\n[red bold]Error:[/red bold] {str(e)}")
        raise typer.Exit(1)

@app.command()
def supported_formats():
    """Show information about supported file formats."""
    table = create_analysis_table("Supported File Formats")
    table.add_column("Format", style="cyan")
    table.add_column("Extensions", style="magenta")
    table.add_column("Description", style="green")
    
    formats = [
        ("FASTA", ".fa, .fasta", "Basic sequence format"),
        ("FASTQ", ".fq, .fastq", "Sequences with quality scores"),
        ("GenBank", ".gb, .gbk, .genbank", "Annotated sequence data"),
        ("EMBL", ".embl", "European sequence format")
    ]
    
    for format_name, extensions, description in formats:
        table.add_row(format_name, extensions, description)
    
    console.print()
    console.print(table)

if __name__ == "__main__":
    app()