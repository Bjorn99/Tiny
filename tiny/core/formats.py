from pathlib import Path
from typing import List, Dict, Union, Optional
from Bio import SeqIO
from Bio import GenBank
import pysam
from dataclasses import dataclass
from Bio.Seq import Seq

@dataclass
class SequenceRecord:
    """Standardize sequence data from different file formats."""
    id: str
    sequence: str
    description: Optional[str] = None
    features: Optional[List[Dict]] = None
    quality_scores: Optional[List[int]] = None
    metadata: Optional[Dict] = None


class FormatHandler:
    
    @staticmethod
    def read_fasta(file_path: Path) -> List[SequenceRecord]:
        """Read sequences from FASTA format."""
        records = []
        for record in SeqIO.parse(str(file_path), "fasta"):
            seq_record = SequenceRecord(
                id=record.id,
                sequence=str(record.seq),
                description=record.description
            )
            records.append(seq_record)
        return records

    @staticmethod
    def read_fastq(file_path: Path) -> List[SequenceRecord]:
        """Read sequences from FASTQ format."""
        records = []
        for record in SeqIO.parse(str(file_path), "fastq"):
            seq_record = SequenceRecord(
                id=record.id,
                sequence=str(record.seq),
                quality_scores=record.letter_annotations["phred_quality"],
                description=record.description
            )
            records.append(seq_record)
        return records

    
    @staticmethod
    def read_genbank(file_path: Path) -> List[SequenceRecord]:
        """Read sequences from GenBank format."""
        records = []
        for record in SeqIO.parse(str(file_path), "genbank"):
            features = []
            for feature in record.features:
                features.append({
                    'type': feature.type,
                    'location': str(feature.location),
                    'qualifiers': feature.qualifiers
                })
            
            seq_record = SequenceRecord(
                id=record.id,
                sequence=str(record.seq),
                description=record.description,
                features=features,
                metadata={
                    'organism': record.annotations.get('organism', ''),
                    'taxonomy': record.annotations.get('taxonomy', []),
                    'references': [ref.title for ref in record.annotations.get('references', [])]
                }
            )
            records.append(seq_record)
        return records

    
    @staticmethod
    def read_embl(file_path: Path) -> List[SequenceRecord]:
        """Read sequences from EMBL format."""
        records = []
        for record in SeqIO.parse(str(file_path), "embl"):
            features = []
            for feature in record.features:
                features.append({
                    'type': feature.type,
                    'location': str(feature.location),
                    'qualifiers': feature.qualifiers
                })
            
            seq_record = SequenceRecord(
                id=record.id,
                sequence=str(record.seq),
                description=record.description,
                features=features,
                metadata={
                    'accessions': record.annotations.get('accessions', []),
                    'taxonomy': record.annotations.get('taxonomy', []),
                    'references': [ref.title for ref in record.annotations.get('references', [])]
                }
            )
            records.append(seq_record)
        return records


    @staticmethod
    def read_sam_bam(file_path: Path) -> List[SequenceRecord]:
        """Read sequences from SAM/BAM format."""
        records = []
        with pysam.AlignmentFile(str(file_path), "r") as sam_file:
            for read in sam_file.fetch():
                if not read.is_unmapped:
                    seq_record = SequenceRecord(
                        id=read.query_name,
                        sequence=read.query_sequence,
                        quality_scores=read.query_qualities,
                        metadata={
                            'mapping_quality': read.mapping_quality,
                            'reference_name': read.reference_name,
                            'reference_start': read.reference_start,
                            'is_reverse': read.is_reverse
                        }
                    )
                    records.append(seq_record)
        return records


    @classmethod
    def read_file(cls, file_path: Path) -> List[SequenceRecord]:
        """
        Read sequences from file based on file extension.
        Supports: FASTA, FASTQ, GenBank, EMBL, SAM, BAM
        """
        suffix = file_path.suffix.lower()
        
        format_handlers = {
            '.fa': cls.read_fasta,
            '.fasta': cls.read_fasta,
            '.fq': cls.read_fastq,
            '.fastq': cls.read_fastq,
            '.gb': cls.read_genbank,
            '.gbk': cls.read_genbank,
            '.genbank': cls.read_genbank,
            '.embl': cls.read_embl,
            '.sam': cls.read_sam_bam,
            '.bam': cls.read_sam_bam
        }
        
        handler = format_handlers.get(suffix)
        if handler is None:
            raise ValueError(f"Unsupported file format: {suffix}")
            
        return handler(file_path)