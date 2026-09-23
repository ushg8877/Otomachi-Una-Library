#!/usr/bin/env python3
"""
Local IO interactor: connect two executables via stdin/stdout pipes.
Does NOT compile — pass already-built binaries.

Protocol:
  - solution stdin  <->  jury stdout
  - solution stdout <->  jury stdin
  - jury stderr = verdict (printed by this tester)

Usage (Linux):
  g++ -O2 -o a a.cpp
  g++ -O2 -o a-jury a-jury.cpp
  python3 interact.py a a-jury
  python3 interact.py ./guess ./guess-jury --in guess.in
  python3 interact.py /path/to/sol /path/to/jury --tl 2
  python3 interact.py a a-jury -q          # quiet, no interaction log
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
import time
from pathlib import Path


def resolve_exe(name: str) -> Path:
    """Resolve an executable path. Bare names look in cwd first (./name, ./name.exe)."""
    path = Path(name)
    candidates: list[Path] = []

    if path.is_absolute():
        candidates.append(path)
    elif name.startswith("./") or name.startswith("../") or os.sep in name or (
        os.altsep and os.altsep in name
    ):
        candidates.append(path)
    else:
        # bare name: prefer cwd (Linux PATH won't find ./a)
        candidates.append(Path.cwd() / path.name)

    # Windows: bare "g" should also try "g.exe"
    extra: list[Path] = []
    for c in candidates:
        if c.suffix.lower() != ".exe":
            extra.append(c.with_suffix(c.suffix + ".exe") if c.suffix else c.with_name(c.name + ".exe"))
    candidates.extend(extra)

    exe = next((c for c in candidates if c.exists()), None)
    if exe is None:
        tried = ", ".join(str(c) for c in candidates)
        raise SystemExit(f"missing executable: {name} (tried: {tried})")
    if not os.access(exe, os.X_OK):
        raise SystemExit(f"not executable: {exe} (chmod +x?)")
    return exe.resolve()


_RESET = "\033[0m"
# < participant (bright yellow), > jury (bright magenta) — high contrast
_COLOR = {"<": "\033[93m", ">": "\033[95m"}
_DIM = "\033[2m"

_print_lock = threading.Lock()


def log_line(marker: str, text: str) -> None:
    color = _COLOR.get(marker, "")
    with _print_lock:
        print(f"{color}{marker}{_RESET} {text}", flush=True)


def log_event(text: str) -> None:
    with _print_lock:
        print(f"{_DIM}* {_RESET}{text}", flush=True)


def pump(src, dst, show: bool, marker: str, done: threading.Event) -> None:
    try:
        while not done.is_set():
            line = src.readline()
            if not line:
                break
            if show:
                text = line.decode("utf-8", errors="replace").rstrip("\r\n")
                log_line(marker, text)
            try:
                dst.write(line)
                dst.flush()
            except (BrokenPipeError, OSError):
                break
    finally:
        try:
            dst.close()
        except Exception:
            pass


def read_stream(src, sink: list[bytes], done: threading.Event) -> None:
    try:
        while not done.is_set():
            chunk = src.read(4096)
            if not chunk:
                break
            sink.append(chunk)
    except Exception:
        pass


def run_interact(
    sol_exe: Path,
    jury_exe: Path,
    jury_args: list[str],
    tl: float | None,
    verbose: bool,
) -> int:
    sol = subprocess.Popen(
        [str(sol_exe)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0,
    )
    jury = subprocess.Popen(
        [str(jury_exe), *jury_args],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0,
    )

    assert sol.stdin and sol.stdout and sol.stderr
    assert jury.stdin and jury.stdout and jury.stderr

    done = threading.Event()
    jury_err: list[bytes] = []
    sol_err: list[bytes] = []

    t1 = threading.Thread(
        target=pump,
        args=(sol.stdout, jury.stdin, verbose, "<", done),
        daemon=True,
    )
    t2 = threading.Thread(
        target=pump,
        args=(jury.stdout, sol.stdin, verbose, ">", done),
        daemon=True,
    )
    t3 = threading.Thread(target=read_stream, args=(jury.stderr, jury_err, done), daemon=True)
    t4 = threading.Thread(target=read_stream, args=(sol.stderr, sol_err, done), daemon=True)

    t0 = time.perf_counter()
    for t in (t1, t2, t3, t4):
        t.start()

    timed_out = False
    sol_marked = False
    jury_marked = False

    def mark_exit(marker: str, code: int | None, noted: bool) -> bool:
        if noted or code is None:
            return noted
        if verbose:
            log_line(marker, f"[exit {code}]")
        return True

    while True:
        sc, jc = sol.poll(), jury.poll()
        if sc is not None:
            sol_marked = mark_exit("<", sc, sol_marked)
        if jc is not None:
            jury_marked = mark_exit(">", jc, jury_marked)

        if sc is not None and jc is not None:
            break
        if jc is not None and sc is None:
            try:
                sol.wait(timeout=0.2)
            except subprocess.TimeoutExpired:
                pass
            if sol.poll() is None:
                if verbose:
                    log_event("sol still running after jury exit; killing sol")
                sol.kill()
            sol_marked = mark_exit("<", sol.poll(), sol_marked)
            break
        if sc is not None and jc is None:
            try:
                jury.wait(timeout=0.2)
            except subprocess.TimeoutExpired:
                pass
            if jury.poll() is None:
                if verbose:
                    log_event("jury still running after sol exit; killing jury")
                jury.kill()
            jury_marked = mark_exit(">", jury.poll(), jury_marked)
            break
        if tl is not None and tl > 0 and time.perf_counter() - t0 > tl:
            timed_out = True
            if verbose:
                log_event(f"TIMEOUT after {tl}s; killing sol and jury")
            sol.kill()
            jury.kill()
            sol_marked = mark_exit("<", sol.poll(), sol_marked)
            jury_marked = mark_exit(">", jury.poll(), jury_marked)
            break
        time.sleep(0.005)

    done.set()
    for t in (t1, t2, t3, t4):
        t.join(timeout=1.0)

    for p in (sol, jury):
        try:
            p.wait(timeout=0.5)
        except subprocess.TimeoutExpired:
            p.kill()

    sol_marked = mark_exit("<", sol.returncode, sol_marked)
    jury_marked = mark_exit(">", jury.returncode, jury_marked)
    elapsed = time.perf_counter() - t0

    sol_err_text = b"".join(sol_err).decode("utf-8", errors="replace").strip()
    jury_err_text = b"".join(jury_err).decode("utf-8", errors="replace").strip()

    if sol_err_text:
        print("---------- solution stderr ----------", flush=True)
        print(sol_err_text, flush=True)
        print("-------------------------------------", flush=True)

    print("---------- jury stderr (result) ----------", flush=True)
    print(jury_err_text if jury_err_text else "(empty)", flush=True)
    print("------------------------------------------", flush=True)

    print(
        f"[status] sol_exit={sol.returncode} jury_exit={jury.returncode} "
        f"time={elapsed:.3f}s" + (" TIMEOUT" if timed_out else ""),
        flush=True,
    )

    if timed_out:
        return 124
    if sol.returncode != 0:
        return 1
    if jury.returncode != 0:
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="IO interactive tester (executables only, no compile)",
        epilog="example: python3 interact.py a a-jury --in test.in",
    )
    ap.add_argument("sol", help="solution executable, e.g. a or ./a")
    ap.add_argument("jury", help="jury executable, e.g. a-jury or ./a-jury")
    ap.add_argument("--in", dest="infile", default=None, help="optional file passed as argv[1] to jury")
    ap.add_argument("--tl", type=float, default=10.0, help="time limit in seconds (default 10; 0=unlimited)")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not print interaction lines")
    ap.add_argument("-C", "--chdir", default=None, help="work directory (default: cwd)")
    args = ap.parse_args()

    if args.chdir:
        os.chdir(args.chdir)

    sol_exe = resolve_exe(args.sol)
    jury_exe = resolve_exe(args.jury)

    jury_args: list[str] = []
    if args.infile:
        infile = Path(args.infile)
        if not infile.exists():
            raise SystemExit(f"missing input file {infile}")
        jury_args.append(str(infile.resolve()))

    print(f"[run] {sol_exe}  <->  {jury_exe} {' '.join(jury_args)}".rstrip(), flush=True)
    return run_interact(sol_exe, jury_exe, jury_args, args.tl, verbose=not args.quiet)


if __name__ == "__main__":
    raise SystemExit(main())
