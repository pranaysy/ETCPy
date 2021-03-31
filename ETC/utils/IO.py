#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains helper functions for reading and writing files.

@author: Pranay S. Yadav
"""

from csv import DictWriter

# Import functions from standard library modules
from pathlib import Path


def populate_files(filepath, suffix="*.txt"):
    if not isinstance(filepath, Path):
        filepath = Path(filepath)

    if filepath.exists() and filepath.is_dir():
        return filepath.rglob(suffix)

    print("Invalid path")
    return None


def read(filepath, delimiter=None):
    if not isinstance(filepath, Path):
        filepath = Path(filepath)
    text = filepath.read_text()

    if delimiter:
        text = "".join(text.split(delimiter))

    return text


def save(out, filename):

    with open(filename, "w") as fileout:
        writer = DictWriter(fileout, fieldnames=out[0].keys(), delimiter=",")
        writer.writeheader()
        writer.writerows(out)
        print(f">> Data successfully stored to disk as {filename}")
