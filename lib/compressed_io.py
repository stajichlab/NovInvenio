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
import io
import subprocess
from pathlib import Path
from typing import IO


class _ZstdReader:
    """Wraps a zstd process's stdout to detect non-zero exit (corrupt/truncated).

    The file may be multi-GB, so we stream it rather than reading into memory.
    The underlying stdout returns an empty stream if zstd fails, but we detect
    that failure by calling proc.wait() on close/exhaustion and raising if
    returncode != 0 (except SIGPIPE, -13/141, which is tolerated when the
    stream is closed early without reading to EOF). This makes a corrupt or
    truncated .zst file fail loudly instead of silently returning zero bytes --
    critical for pipelines where a meaningless empty result is worse than a
    crash (e.g. column ordering by locus position, where alphabetical fallback
    is wrong but plausible-looking).

    Supports context manager, iteration, read(), readline(), readlines(), and
    csv.DictReader use.
    """

    def __init__(self, proc: subprocess.Popen) -> None:
        self._proc = proc
        self._stdout = proc.stdout
        self._stderr = proc.stderr
        self._text_wrapper = None
        # Give zstd a moment to fail on corrupt input
        # (it fails fast, so poll() will get the error immediately)
        # Only do this check early if the process has already exited
        if proc.poll() is not None and proc.returncode != 0:
            stderr = self._get_stderr()
            raise RuntimeError(
                f"zstd decompression failed with exit code {proc.returncode}"
                f"{f': {stderr}' if stderr else ''}"
            )

    def _get_stderr(self) -> str:
        """Read stderr if available."""
        if self._stderr:
            try:
                return self._stderr.read().decode("utf-8", errors="replace")
            except Exception:
                return ""
        return ""

    def _ensure_wrapper(self):
        """Lazily create TextIOWrapper on first use, after checking return code."""
        if self._text_wrapper is None:
            try:
                self._text_wrapper = io.TextIOWrapper(self._stdout, encoding="utf-8")
            except UnicodeDecodeError:
                # Corrupt data - check if zstd failed
                self._proc.wait()
                if self._proc.returncode != 0:
                    stderr = self._get_stderr()
                    raise RuntimeError(
                        f"zstd decompression failed with exit code {self._proc.returncode}"
                        f"{f': {stderr}' if stderr else ''}"
                    ) from None
                raise
        return self._text_wrapper

    def _safe_read(self, func, *args):
        """Call a read function, raising RuntimeError if zstd failed."""
        try:
            return func(*args)
        except UnicodeDecodeError as e:
            # Text decoding failed - likely corrupt zstd output
            # Wait for the process to complete
            returncode = self._proc.wait()
            stderr = self._get_stderr()
            if returncode != 0:
                # Process exited with error - report that
                raise RuntimeError(
                    f"zstd decompression failed with exit code {returncode}"
                    f"{f': {stderr}' if stderr else ''}"
                ) from None
            # Even if zstd exited 0, corrupt output that can't be decoded is an error
            raise RuntimeError(
                "zstd decompression produced invalid UTF-8 output (corrupt input?)"
            ) from e

    def read(self, size: int = -1) -> str:
        wrapper = self._ensure_wrapper()
        return self._safe_read(wrapper.read, size)

    def readline(self) -> str:
        wrapper = self._ensure_wrapper()
        return self._safe_read(wrapper.readline)

    def readlines(self) -> list[str]:
        wrapper = self._ensure_wrapper()
        return self._safe_read(wrapper.readlines)

    def __iter__(self):
        return self

    def __next__(self) -> str:
        wrapper = self._ensure_wrapper()
        try:
            return wrapper.__next__()
        except UnicodeDecodeError:
            returncode = self._proc.poll()
            if returncode is not None and returncode != 0:
                stderr = self._get_stderr()
                raise RuntimeError(
                    f"zstd decompression failed with exit code {returncode}"
                    f"{f': {stderr}' if stderr else ''}"
                ) from None
            raise

    def __enter__(self) -> "_ZstdReader":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Close stdout and check zstd's exit code. If it failed, raise.

        SIGPIPE (-13 or exit code 141) is tolerated because it occurs when
        the stream is closed early without reading to EOF, which is a normal
        use case (e.g., breaking from a loop). Any other non-zero exit is
        treated as a decompression error.
        """
        if self._text_wrapper:
            self._text_wrapper.close()
        else:
            self._stdout.close()
        returncode = self._proc.wait()
        # Tolerate SIGPIPE (return code -13 appears as 141 after wait())
        if returncode != 0 and returncode not in (-13, 141):
            stderr = self._get_stderr()
            raise RuntimeError(
                f"zstd decompression failed with exit code {returncode}"
                f"{f': {stderr}' if stderr else ''}"
            )


def open_maybe_compressed(path: str | Path) -> IO[str]:
    """Open `path` for text reading, transparently decompressing based on
    its extension (.gz, .zst) or returning a plain text handle otherwise.

    All formats check for file existence first, so missing files raise
    FileNotFoundError consistently (not silent empty streams).

    The .zst case shells out to `zstd -dc` (no `-f` flag) and wraps stdout
    in a reader that detects non-zero exit on close/exhaustion. This catches
    corrupt or truncated .zst files that would otherwise silently return empty
    streams. Symlinks (Nextflow staging, PR #113) are resolved to their real
    path in Python BEFORE passing to zstd, so zstd never sees a symlink and
    doesn't need `-f`. The `-f` flag was tried initially but rejected: it
    silently passes non-zstd files through with exit 0, defeating the purpose
    of this function (fail loudly on invalid input).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"no such file: {path}")

    path_str = str(path)
    if path_str.endswith(".gz"):
        return gzip.open(path_str, "rt")
    if path_str.endswith(".zst"):
        # Resolve symlinks so zstd never sees them and doesn't need -f
        resolved_path = str(p.resolve())
        proc = subprocess.Popen(
            ["zstd", "-dc", resolved_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return _ZstdReader(proc)
    return open(path_str)


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
