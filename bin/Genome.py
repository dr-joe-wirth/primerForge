from __future__ import annotations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Generator


class InvalidFileFormat(Exception):
    pass


class Genome:
    def __init__(self, fn:Path, preload:bool=False) -> None:
        self.fn:Path
        self.length:int
        self.name:str
        self._format:str|None
        self._seqs:list[SeqRecord]
        self._num_seqs:int

        self.fn = fn
        self.length = 0
        self.name = fn.name
        self._num_seqs = 0
        self._seqs = list()
        self._introspect()

        if preload:
            self._load()

    def __iter__(self) -> Generator[SeqRecord,None,None]:
        if not self._seqs:
            for seq in self._stream():
                yield seq

        for seq in self._seqs:
            yield seq

    def __len__(self) -> int:
        return self._num_seqs

    def __eq__(self, other:Genome) -> bool:
        return self.fn.absolute() == other.fn.absolute()

    def __ne__(self, other:Genome) -> bool:
        return not self == other

    def __gt__(self, other:Genome) -> bool:
        return self.length > other.length

    def __hash__(self) -> int:
        return hash(self.fn.absolute())

    def _load(self) -> None:
        if self._format:
            self._seqs = list(self._stream())

    def _stream(self) -> Generator[SeqRecord,None,None]:
        if self._format:
            with open(self.fn, 'r') as h:
                for rec in SeqIO.parse(h, self._format):
                    yield rec

    def _clear(self) -> None:
        self._seqs = list()

    def _sniff_format(self) -> None:
        with open(self.fn, 'r') as h:
            line = h.readline()

            if line.startswith('>'):
                self._format = 'fasta'

            elif line.startswith('LOCUS'):
                self._format = 'genbank'

            else:
                self._format = None

    def _introspect(self) -> None:
        self._sniff_format()

        if self._format:
            for seq in self:
                self.length += len(seq)
                self._num_seqs += 1
