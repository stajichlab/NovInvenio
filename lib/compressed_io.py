"""Transparent gzip/zstd input support for this study's `bin/` scripts --
per the general storage-compression convention (large text intermediates
should default to compressed, and scripts that read them should accept
compressed input directly rather than requiring a manual decompression
step first).

`.gz` via the stdlib `gzip` module; `.zst` by shelling out to the `zstd`
CLI (no extra Python dependency needed -- `zstd` is already a pixi
dependency for this workspace). Plain, uncompressed files still work
unchanged.
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import IO


def open_maybe_compressed(path: str | Path) -> IO[str]:
    """Open `path` for text reading, transparently decompressing based on
    its extension (.gz, .zst) or returning a plain text handle otherwise.

    The .zst case shells out to `zstd -dc` and returns its stdout as a
    text stream -- there is no line-iteration difference from a normal
    file handle, but the subprocess is not explicitly waited on here;
    for a short-lived script that reads the whole stream and exits, the
    OS reaps it on process exit, which is an acceptable simplification
    for how this helper is actually used (one big streaming read, not a
    long-lived server)."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        proc = subprocess.Popen(
            ["zstd", "-dc", path], stdout=subprocess.PIPE, text=True,
        )
        return proc.stdout
    return open(path)


def open_maybe_compressed_write(path: str | Path) -> IO[str]:
    """Open `path` for text writing, transparently compressing based on its
    extension (.gz, .zst) or writing plain text otherwise -- the write-side
    counterpart to `open_maybe_compressed`. Callers pick compression by
    naming the output path with a `.gz`/`.zst` suffix; nothing else about
    the call site changes (`with open_maybe_compressed_write(path) as fh:
    fh.write(...)` reads identically to plain `open`).

    The .zst case pipes through `zstd -q -f` (single-threaded -- this
    writes small-to-medium per-call outputs like a presence matrix, not a
    multi-GB tblastn stream; a multi-GB stream should be compressed via a
    dedicated `zstd -T0` pipe at the shell level instead, as
    run_rescue_pass_tblastn_chunked.sh already does for its own tblastn
    output). The subprocess IS waited on and its exit code checked here
    (unlike the read side): a write that silently produced a truncated or
    empty compressed file because the subprocess died would be far worse
    than a slow decompression, since there would be no second read to
    surface the problem.
    """
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    if path.endswith(".zst"):
        return _ZstdWriter(path)
    return open(path, "w")


class _ZstdWriter:
    """Text-mode write handle that pipes everything written to it through
    `zstd -q -f -o <path>` on close. A context manager (not a plain file
    object) because the compressing subprocess's stdin must be closed
    before its exit code can be checked -- `__exit__` does both, in order,
    so a failed compression raises before the caller can mistake a
    truncated/missing .zst file for a successful write.
    """

    def __init__(self, path: str) -> None:
        self._proc = subprocess.Popen(
            ["zstd", "-q", "-f", "-o", path], stdin=subprocess.PIPE, text=True,
        )

    def write(self, data: str) -> int:
        return self._proc.stdin.write(data)

    def __enter__(self) -> "_ZstdWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._proc.stdin.close()
        returncode = self._proc.wait()
        if exc_type is None and returncode != 0:
            raise subprocess.CalledProcessError(returncode, ["zstd", "-q", "-f", "-o"])
