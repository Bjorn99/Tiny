from pathlib import Path
from typing import List, Dict, Any
import json
from Bio import SeqIO
from tiny.core.sequence import DNASequence

def load_fasta(file_path: Path) -> List[str]:
    """Load sequences from a FASTA file."""
    sequences = []
    for record in SeqIO.parse(str(file_path), "fasta"):
        sequences.append(str(record.seq))
    return sequences

def save_results(results: List[Any], output_path: Path) -> None:
    """Save analysis results to a file."""
    output_data = []
    for result in results:
        result_dict = {
            "sequence": str(result.sequence),
            "length": result.length,
            "gc_content": result.gc_content,
            "molecular_weight": result.molecular_weight,
            "base_composition": result.base_composition
        }
        output_data.append(result_dict)
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)

def format_sequence(sequence: str, width: int = 60) -> str:
    """Format a sequence with line breaks for better readability."""
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))