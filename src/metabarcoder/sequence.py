from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator


DNA_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict")
    return path.open("rt", encoding="ascii", errors="strict")


def iter_fastq(paths: list[Path]) -> Iterator[tuple[str, str, str]]:
    for path in paths:
        with open_text(path) as handle:
            line_number = 0
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip("\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\r\n")
                line_number += 4
                if not sequence or not plus or not quality:
                    raise ValueError(f"Neúplný FASTQ záznam v {path} pri riadku {line_number}")
                if not header.startswith("@") or not plus.startswith("+"):
                    raise ValueError(f"Neplatný FASTQ formát v {path} pri riadku {line_number - 3}")
                if len(sequence) != len(quality):
                    raise ValueError(f"Sekvencia a kvalita majú inú dĺžku v {path}")
                yield header[1:].strip(), sequence.upper(), quality


def mean_phred(quality: str) -> float:
    return sum(ord(char) - 33 for char in quality) / len(quality) if quality else 0.0


def levenshtein(a: str, b: str, max_distance: int | None = None) -> int:
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for i, char_b in enumerate(b, 1):
        current = [i]
        row_min = i
        for j, char_a in enumerate(a, 1):
            value = min(
                current[j - 1] + 1,
                previous[j] + 1,
                previous[j - 1] + (char_a != char_b),
            )
            current.append(value)
            row_min = min(row_min, value)
        previous = current
        if max_distance is not None and row_min > max_distance:
            return max_distance + 1
    return previous[-1]


def approximate_substring(sequence: str, pattern: str, max_errors: int) -> tuple[int, int, int] | None:
    """Return (start, end, edits) for the best approximate pattern occurrence.

    This is a semi-global edit-distance alignment: the complete primer must
    align, while unaligned sequence before and after it has no penalty. It
    tolerates substitutions and short indels without testing every possible
    window separately.
    """
    sequence = sequence.upper()
    pattern = pattern.upper()
    if not pattern or not sequence:
        return None
    exact_start = sequence.find(pattern)
    if exact_start >= 0:
        return exact_start, exact_start + len(pattern), 0
    width = len(sequence) + 1
    previous_cost = [0] * width
    previous_start = list(range(width))
    for row, pattern_base in enumerate(pattern, 1):
        current_cost = [row] + [0] * len(sequence)
        current_start = [0] + [0] * len(sequence)
        for column, sequence_base in enumerate(sequence, 1):
            candidates = (
                (previous_cost[column] + 1, previous_start[column]),
                (current_cost[column - 1] + 1, current_start[column - 1]),
                (previous_cost[column - 1] + (pattern_base != sequence_base),
                 previous_start[column - 1]),
            )
            cost, start = min(candidates, key=lambda item: (item[0], -item[1]))
            current_cost[column] = cost
            current_start[column] = start
        previous_cost, current_cost = current_cost, previous_cost
        previous_start, current_start = current_start, previous_start
    end = min(range(1, width), key=lambda column: previous_cost[column])
    distance = previous_cost[end]
    if distance > max_errors:
        return None
    return previous_start[end], end, distance


def find_amplicon(sequence: str, forward: str, reverse_primer: str, max_errors: int,
                   min_insert: int, max_insert: int) -> tuple[str, bool] | None:
    reverse_binding = reverse_complement(reverse_primer)
    for oriented, was_reversed in ((sequence, False), (reverse_complement(sequence), True)):
        left = approximate_substring(oriented, forward, max_errors)
        if left is None:
            continue
        search_start = left[1] + min_insert
        search_end = min(len(oriented), left[1] + max_insert + len(reverse_binding) + max_errors)
        if search_start >= len(oriented):
            continue
        right_local = approximate_substring(oriented[search_start:search_end], reverse_binding, max_errors)
        if right_local is None:
            continue
        right_start = search_start + right_local[0]
        insert = oriented[left[1]:right_start]
        if min_insert <= len(insert) <= max_insert:
            return insert, was_reversed
    return None


def discover_fastq(inputs: list[Path]) -> list[Path]:
    discovered: list[Path] = []
    for item in inputs:
        if item.is_dir():
            discovered.extend(path for path in item.rglob("*") if path.is_file() and
                              (path.name.lower().endswith(".fastq") or path.name.lower().endswith(".fastq.gz")))
        elif item.is_file() and (item.name.lower().endswith(".fastq") or item.name.lower().endswith(".fastq.gz")):
            discovered.append(item)
        else:
            raise ValueError(f"Nepodporovaný alebo chýbajúci vstup: {item}")
    unique = sorted(set(path.resolve() for path in discovered))
    if not unique:
        raise ValueError("Nenašiel sa žiadny FASTQ alebo FASTQ.gz súbor")
    return unique

