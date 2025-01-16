import typer
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from typing import List, Optional
from pathlib import Path
from enum import Enum

from tiny.core.sequence import DNASequence
from tiny.core.analysis import analyze_sequences, compare_sequences
from tiny.core.advanced import AdvancedAnalysis
from tiny.core.utils import load_fasta, save_results

app = typer.Typer()
console = Console()

class AlignmentMode(str, Enum):
    global_align = "global"
    local = "local"
    semi_global = "semi-global"

@app.command()
def analyze(
    sequences: List[str] = typer.Argument(None),
    fasta_file: Optional[Path] = typer.Option(None, "--fasta", "-f", help="Input FASTA file"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output file for results"),
):
    """Analyze DNA sequences for GC content and other properties."""
    if fasta_file:
        sequences = load_fasta(fasta_file)
    elif not sequences:
        console.print("[red]Error: Please provide sequences or a FASTA file[/red]")
        raise typer.Exit(1)

    results = analyze_sequences([DNASequence(seq) for seq in sequences])
    
    # Create a rich table for display
    table = Table(title="Sequence Analysis Results")
    table.add_column("Sequence", style="cyan")
    table.add_column("Length", style="magenta")
    table.add_column("GC Content", style="green")
    table.add_column("Molecular Weight", style="yellow")

    for result in results:
        table.add_row(
            str(result.sequence)[:30] + "..." if len(str(result.sequence)) > 30 else str(result.sequence),
            str(result.length),
            f"{result.gc_content:.2f}%",
            f"{result.molecular_weight:.2f}"
        )

    console.print(table)
    
    if output:
        save_results(results, output)

@app.command()
def compare(
    sequence1: str = typer.Argument(..., help="First DNA sequence"),
    sequence2: str = typer.Argument(..., help="Second DNA sequence"),
):
    """Compare two DNA sequences for mutations and similarities."""
    comparison = compare_sequences(DNASequence(sequence1), DNASequence(sequence2))
    
    console.print("\n[bold]Sequence Comparison Results[/bold]")
    console.print(f"Identity: [green]{comparison.identity:.2f}%[/green]")
    console.print(f"Alignment Score: [blue]{comparison.alignment_score}[/blue]")
    
    if comparison.mutations:
        console.print("\n[bold red]Mutations:[/bold red]")
        for mut in comparison.mutations:
            console.print(f"Position {mut.position}: {mut.original} → {mut.mutated}")

@app.command()
def align(
    sequence1: str = typer.Argument(..., help="First DNA sequence"),
    sequence2: str = typer.Argument(..., help="Second DNA sequence"),
    mode: AlignmentMode = typer.Option(
        AlignmentMode.global_align,
        help="Alignment mode"
    ),
):
    """sequence alignment with detailed visualization."""
    try:
        result = AdvancedAnalysis.align_sequences(sequence1, sequence2, mode)
        
        # Rich text visualization
        console.print("\n[bold]Sequence Alignment Results[/bold]")
        console.print(f"Mode: [cyan]{mode}[/cyan]")
        console.print(f"Alignment Score: [green]{result.alignment_score:.2f}[/green]")
        console.print(f"Sequence Identity: [blue]{result.identity:.2f}%[/blue]")
        console.print(f"Gaps: [red]{result.gaps}[/red]")
        
        # Visualize alignment
        console.print("\n[bold]Alignment Visualization:[/bold]")
        console.print(f"Sequence 1: [cyan]{result.aligned_seq1}[/cyan]")
        console.print(f"Sequence 2: [cyan]{result.aligned_seq2}[/cyan]")
        
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")

@app.command()
def find_motifs(
    sequences: List[str] = typer.Argument(None),
    fasta_file: Optional[Path] = typer.Option(None, "--fasta", "-f"),
    motif_length: int = typer.Option(6, "--length", "-l"),
    min_frequency: int = typer.Option(2, "--min-freq", "-m"),
):
    """Find common motifs in DNA sequences."""
    if fasta_file:
        sequences = load_fasta(fasta_file)
    elif not sequences:
        console.print("[red]Error: Please provide sequences or a FASTA file[/red]")
        raise typer.Exit(1)

    try:
        motifs = AdvancedAnalysis.find_motifs(sequences, motif_length, min_frequency)
        
        # Results table
        table = Table(title="Motif Analysis Results")
        table.add_column("Motif", style="cyan")
        table.add_column("Frequency", style="magenta")
        table.add_column("Consensus Score", style="green")
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
        console.print(f"[red]Error: {str(e)}[/red]")

@app.command()
def find_regulatory(
    sequence: str = typer.Argument(..., help="DNA sequence to analyze"),
):
    """Find potential regulatory elements in a DNA sequence."""
    try:
        elements = AdvancedAnalysis.find_regulatory_elements(sequence)
        
        console.print("\n[bold]Regulatory Elements Analysis[/bold]")
        
        for element_type, findings in elements.items():
            if findings:
                panel_content = "\n".join(
                    f"Position {pos}: {seq}" for pos, seq in findings
                )
                console.print(Panel(
                    panel_content,
                    title=f"[cyan]{element_type.replace('_', ' ').title()}[/cyan]",
                    expand=False
                ))
            else:
                console.print(f"\nNo {element_type.replace('_', ' ').title()} found.")

    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")

if __name__ == "__main__":
    app()