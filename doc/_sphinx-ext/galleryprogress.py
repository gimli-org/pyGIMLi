"""Condensed progress report for the documentation gallery.

The examples print their modelling and inversion history from the C++ core of
pyGIMLi. That output bypasses the Python level capturing of sphinx-gallery and
floods the build log, which makes it impossible to see how far the gallery is.

Here every example gets its slot of a known list. Its stdout is redirected
into a log file for the duration of the example, and the state of all slots is
kept in ``_gallery_status`` next to the gallery. In a terminal a finished
example is written as one permanent line and the running ones are kept in a
small table at the end of the terminal::

    [ 12/57] examples/3_ert/plot_02_ert_field_data.py ... [ok]      7.9 s
    [ 13/57] examples/2_seismics/plot_04_koenigsee.py ... [ok]     21.4 s
    [ 18/57] examples/3_ert/plot_07_ert_3dsurface.py ... [building]
    [ 19/57] examples/3_ert/plot_08_ert_crosshole3d.py ... [building]
    13/57 done, 2 running

While the table is active the file descriptors 1 and 2 of the build process
are routed through a pipe, so foreign output (Sphinx warnings, logging, C++
messages) is written above the table instead of into it. Without a terminal
(e.g. ``--target doc > log``) one line per state change is written instead.
The console output of an example is kept in ``_gallery_logs/<example>.log``.

The gallery has to be configured with ``reset_modules_order = "both"`` so that
the hook is called before and after each example. Set ``GIMLI_DOC_PROGRESS=0``
to get the unfiltered output back.
"""

import fcntl
import json
import os
import re
import select
import shutil
import sys
import threading
import time
import warnings
from pathlib import Path

__all__ = [
    "exampleTimes",
    "finishProgress",
    "flushProgress",
    "markFailed",
    "markWarmed",
    "reportProgress",
    "resetProgress",
    "setup",
]

_enabled = os.environ.get("GIMLI_DOC_PROGRESS", "1") != "0"
# The table is redrawn in place, so its block has to stay at the end of the
# terminal. That holds because all output of the build process is routed
# through a pipe; set GIMLI_DOC_LIVE=0 to get plain lines instead.
_live = os.environ.get("GIMLI_DOC_LIVE", "1") != "0"
_examples = None
_current = None
_printer = None
_warnings = []

_NAME_WIDTH = 70
# The table only holds the running examples, so the block stays small.
_MAX_ROWS = 24
_GREEN = "\033[32m"
_RED = "\033[31m"
_YELLOW = "\033[33m"
_RESET = "\033[0m"


def _galleryDir(galleryConf):
    """Return the output directory of the first gallery."""
    galleryDirs = galleryConf.get("gallery_dirs") or ["_gallery"]
    srcDir = Path(str(galleryConf.get("src_dir") or "."))
    return srcDir / str(galleryDirs[0])


def _statusDir(galleryConf):
    """Return the directory that holds one state file per example slot."""
    return _galleryDir(galleryConf) / "_gallery_status"


def _liveMarker(galleryConf):
    """Return the file that marks an active table in the Sphinx process."""
    return _statusDir(galleryConf) / "_table"


def _widthFile(galleryConf):
    """Return the file that holds the name width of the running table."""
    return _statusDir(galleryConf) / "_width"


def _nameWidth(galleryConf):
    """Return the name width that fits into the terminal of the build.

    The examples run in worker processes that cannot ask the terminal of the
    build, so the process that draws the table tells them the width. A line
    wider than the terminal would wrap and break the redrawing of the table.
    """
    try:
        width = int(_widthFile(galleryConf).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return _NAME_WIDTH
    return max(20, min(_NAME_WIDTH, width))


def _isLive(galleryConf):
    """Return whether the Sphinx process draws the table for this gallery."""
    return _enabled and _liveMarker(galleryConf).exists()


def _warmedMarker(galleryConf):
    """Return the file that marks all examples as executed already."""
    return _statusDir(galleryConf) / "_warmed"


def _isWarmed(galleryConf):
    """Return whether the examples were executed in the global pool."""
    return _warmedMarker(galleryConf).exists()


def resetProgress(galleryConf):
    """Forget the states of a previous run.

    Only the marker of a running table is kept: the warm-up calls this after
    the table was started, and a deleted marker would make the workers believe
    they are not live and print their own lines. The ``_warmed`` marker has to
    go, otherwise the hooks of this build believe the examples are done and the
    example output ends up on the terminal again.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration.
    """
    if not _enabled:
        return
    directory = _statusDir(galleryConf)
    if not directory.is_dir():
        return
    keep = _liveMarker(galleryConf).name
    for path in directory.iterdir():
        if path.name == keep:
            continue
        if path.name.isdigit() or path.name.startswith("_"):
            try:
                path.unlink()
            except OSError:
                continue


def markWarmed(galleryConf):
    """Stop the progress hooks, the examples are executed already.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration.
    """
    if not _enabled:
        return
    directory = _statusDir(galleryConf)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        _warmedMarker(galleryConf).write_text("warmed", encoding="utf-8")
    except OSError:
        pass


def _orderedExamples(galleryConf):
    """Return the example names in the order sphinx-gallery runs them.

    The names carry the name of the gallery directory, so examples and
    tutorials of the same subsection can be told apart, e.g.
    ``examples/3_ert/plot_02_ert_field_data.py``.
    """
    includeRe = re.compile(str(galleryConf.get("filename_pattern", "/plot_")))
    ignoreRe = re.compile(str(galleryConf.get("ignore_pattern", r"\.\.\.")))
    srcDir = Path(str(galleryConf.get("src_dir") or "."))
    names = []
    for exampleDir in galleryConf.get("examples_dirs") or []:
        root = srcDir / str(exampleDir)
        label = Path(str(exampleDir)).name
        for path in sorted(root.rglob("*.py")):
            relative = path.relative_to(root).as_posix()
            if not includeRe.search(f"/{relative}"):
                continue
            if ignoreRe.search(f"/{path.as_posix()}"):
                continue
            names.append(f"{label}/{relative}")
    return names


def _slot(galleryConf, fname):
    """Return the table slot (1 based) of *fname*, 0 if it is unknown."""
    global _examples

    if _examples is None:
        names = _orderedExamples(galleryConf)
        slots = {Path(name).name: i + 1 for i, name in enumerate(names)}
        _examples = (names, slots)

    _, slots = _examples
    return slots.get(Path(str(fname)).name, 0)


def _displayName(index, fname):
    """Return the prefixed example name used in the output."""
    if index and _examples:
        return _examples[0][index - 1]
    return Path(str(fname)).name


def _shorten(name, width):
    """Cut an example name to *width*, the file name is the useful part."""
    if len(name) <= width:
        return name
    if width <= 3:
        return name[-width:]
    return "..." + name[-(width - 3) :]


def flushProgress():
    """Report the examples that finished since the last refresh.

    The table is refreshed by a polling thread, so sphinx can be faster than
    the next refresh. Without this the line of the last example would appear
    after the output that sphinx writes next.
    """
    printer = _printer
    if printer is None or not printer.is_alive():
        return
    try:
        printer.render()
    except OSError:
        pass


def finishProgress():
    """Remove the block of the running examples and report the warnings.

    Called when all examples have been executed, the rest of the build should
    not have a status line for the gallery in the terminal anymore.
    """
    printer = _printer
    if printer is not None:
        try:
            printer.finish()
        except OSError:
            pass
    _reportWarnings()


def _reportWarnings():
    """Print the warnings of the build that were collected meanwhile."""
    if not _warnings:
        return
    count = len(_warnings)
    print(
        f"gallery: {count} warning{'' if count == 1 else 's'} "
        f"while the examples ran",
        flush=True,
    )
    for name, message, source, lineno in _warnings:
        print(f"  {name}: {message} ({source}:{lineno})", flush=True)
    _warnings.clear()


def _setStatus(galleryConf, index, state):
    """Store the state of one example, one file per slot needs no locking."""
    if not index:
        return
    directory = _statusDir(galleryConf)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        (directory / str(index)).write_text(state, encoding="utf-8")
    except OSError:
        pass


def _timeFile(galleryConf):
    """Return the file that holds the runtime of the examples.

    It lives next to the gallery directories, so that a ``clean-gallery``
    does not remove the runtimes the scheduling of the next run is based on.
    """
    srcDir = Path(str(galleryConf.get("src_dir") or "."))
    return srcDir / "_gallery_times.json"


def _storeTime(galleryConf, fname, duration):
    """Store the runtime of one example for the scheduling of the next run."""
    path = _timeFile(galleryConf)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a+", encoding="utf-8") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                handle.seek(0)
                try:
                    times = json.loads(handle.read() or "{}")
                except ValueError:
                    times = {}
                times[Path(str(fname)).name] = round(duration, 3)
                handle.seek(0)
                handle.truncate()
                json.dump(times, handle)
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    except OSError:
        pass


def exampleTimes(galleryConf):
    """Return the cached runtime per example name.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration.

    Returns
    -------
    dict
        Mapping of the example file name to its runtime in seconds.
    """
    try:
        return json.loads(_timeFile(galleryConf).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def _readStatus(galleryConf):
    """Return the states of all examples that have been started."""
    states = {}
    directory = _statusDir(galleryConf)
    if not directory.is_dir():
        return states
    for path in directory.iterdir():
        try:
            states[int(path.name)] = path.read_text(encoding="utf-8").strip()
        except (OSError, ValueError):
            continue
    return states


class _TablePrinter(threading.Thread):
    """Keep the running examples at the end of the terminal.

    Finished examples are written as plain lines, the running ones as a small
    block that is erased and redrawn on every change. To keep the block at the
    end of the terminal, the file descriptors 1 and 2 of the build process are
    routed through a pipe: everything else that is printed (Sphinx warnings,
    logging, C++ messages) is forwarded above the block instead of into it.
    """

    def __init__(self, galleryConf, names):
        super().__init__(daemon=True)
        self.galleryConf = galleryConf
        self.names = names
        self.total = len(names)
        self.lock = threading.RLock()
        self.stopEvent = threading.Event()
        self.terminalFd = os.dup(1)
        self.readFd = None
        self.printedLines = 0
        self.pending = ""
        self.lastStates = None
        self.blockVisible = True
        self.columns = 80
        self.room = _MAX_ROWS
        self.nameWidth = _NAME_WIDTH
        self.updateSize()

    def updateSize(self):
        """Adapt the layout to the size of the terminal."""
        try:
            size = os.get_terminal_size(self.terminalFd)
        except OSError:
            size = shutil.get_terminal_size()
        self.columns = size.columns
        self.room = max(3, min(_MAX_ROWS, size.lines - 2))
        # Rows must not wrap, otherwise the erase would move the cursor to the
        # wrong line and the block would be written a second time. A result
        # line needs nameWidth + 30 columns, two columns are kept free.
        self.nameWidth = max(20, min(_NAME_WIDTH, size.columns - 32))
        self.storeWidth()

    def storeWidth(self):
        """Tell the worker processes how wide their output lines may be."""
        path = _widthFile(self.galleryConf)
        try:
            current = path.read_text(encoding="utf-8").strip()
        except OSError:
            current = ""
        if current == str(self.nameWidth):
            return
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(str(self.nameWidth), encoding="utf-8")
        except OSError:
            pass

    def capture(self):
        """Route fd 1 and fd 2 through a pipe that feeds :meth:`writeText`."""
        readFd, writeFd = os.pipe()
        os.dup2(writeFd, 1)
        os.dup2(writeFd, 2)
        os.close(writeFd)
        self.readFd = readFd

    def release(self):
        """Send the descriptors back to the terminal and drain the pipe."""
        os.dup2(self.terminalFd, 1)
        os.dup2(self.terminalFd, 2)
        self.stopEvent.set()
        self.join(timeout=10)

    def run(self):
        """Forward the captured output and refresh the block."""
        try:
            self.forward()
        finally:
            try:
                # Report the states of examples that finished right at the end.
                self.render()
            except OSError:
                pass
            try:
                with self.lock:
                    self.erase()
                    if self.pending:
                        self.writeLine(self.pending)
                        self.pending = ""
            except OSError:
                pass
            # Never leave the process with a dead pipe on fd 1.
            for number in (1, 2):
                try:
                    os.dup2(self.terminalFd, number)
                except OSError:
                    pass
            try:
                os.close(self.readFd)
            except (OSError, TypeError):
                pass
            self.readFd = None

    def forward(self):
        """Read the pipe until it is closed and refresh the block."""
        while self.readFd is not None:
            if select.select([self.readFd], [], [], 0.2)[0]:
                try:
                    chunk = os.read(self.readFd, 4096)
                except OSError:
                    return
                if not chunk:
                    return
                self.writeText(chunk.decode("utf-8", "replace"))
                continue
            if self.stopEvent.is_set():
                return
            try:
                self.render()
            except OSError:
                return

    def writeText(self, text):
        """Write foreign output above the block."""
        with self.lock:
            self.pending += text
            if "\n" not in self.pending:
                return
            lines = self.pending.split("\n")
            self.pending = lines.pop()
            self.erase()
            for line in lines:
                self.writeLine(line)
            self.draw(self.snapshot())

    def finish(self):
        """Remove the block for good, the examples have been executed."""
        with self.lock:
            self.erase()
            self.blockVisible = False

    def snapshot(self):
        """Return the state of every example slot."""
        states = _readStatus(self.galleryConf)
        return tuple(states.get(i, "") for i in range(1, self.total + 1))

    def render(self):
        """Redraw the block and write the lines of finished examples."""
        snapshot = self.snapshot()
        if snapshot == self.lastStates:
            return
        previous = self.lastStates
        self.lastStates = snapshot
        self.updateSize()
        with self.lock:
            self.erase()
            for index, state in enumerate(snapshot, start=1):
                if not state.startswith(("ok", "fail")):
                    continue
                if previous is not None and previous[index - 1] == state:
                    continue
                self.writeLine(self.finishedLine(index, state))
            self.draw(snapshot)

    def finishedLine(self, index, state):
        """Return the permanent line of a finished example."""
        width = self.nameWidth
        name = _shorten(self.names[index - 1], width)
        head = f"[{index:3d}/{self.total}] {name:<{width}} ..."
        if state.startswith("fail"):
            return f"{head} [{_RED}fail{_RESET}]"
        return f"{head} [{_GREEN}ok{_RESET}]{state[2:]}"

    def writeLine(self, line):
        """Append one permanent line, the caller holds the lock."""
        os.write(self.terminalFd, f"{line}\r\n".encode())

    def erase(self):
        """Remove the block from the terminal."""
        if not self.printedLines:
            return
        # The block is written without a trailing newline, so the cursor sits
        # on its last line. Only CUU (A) and ED (J) are used, CPL (F) is not
        # supported by every terminal emulator.
        text = "\r"
        if self.printedLines > 1:
            text = f"\033[{self.printedLines - 1}A\r"
        os.write(self.terminalFd, f"{text}\033[J".encode())
        self.printedLines = 0

    def draw(self, snapshot):
        """Write the block below the current output."""
        if not self.blockVisible:
            return
        lines = self.body(snapshot)[-self.room :]
        # Every row starts in the first column explicitly, a line feed alone
        # only moves the cursor down.
        text = "\r\033[J" + "\r\n".join(lines)
        os.write(self.terminalFd, text.encode())
        self.printedLines = len(lines)

    def body(self, snapshot):
        """Return the block lines: the running examples and a footer."""
        done = sum(1 for state in snapshot if state.startswith("ok"))
        failed = sum(1 for state in snapshot if state.startswith("fail"))
        running = [
            i
            for i, state in enumerate(snapshot, start=1)
            if state.startswith("building")
        ]
        rows = [
            self.row(i, self.names[i - 1], snapshot[i - 1])
            for i in running[-self.room + 1 :]
        ]
        footer = f"{done}/{self.total} done, {len(running)} running"
        if failed:
            footer += f", {_RED}{failed} failed{_RESET}"
        return rows + [footer]

    def row(self, index, name, state):
        """Return the line of one example."""
        width = self.nameWidth
        shown = _shorten(name, width)
        head = f"[{index:3d}/{self.total}] {shown:<{width}} ..."
        if state.startswith("ok"):
            return f"{head} [{_GREEN}ok{_RESET}]{state[2:]}"
        if state.startswith("fail"):
            return f"{head} [{_RED}fail{_RESET}]"
        if state.startswith("building"):
            return f"{head} [{_YELLOW}building{_RESET}]"
        return f"{head} [pending]"


def _statusText(state):
    """Return the coloured status token of one example."""
    if state == "ok":
        return f"[{_GREEN}ok{_RESET}]"
    if state == "fail":
        return f"[{_RED}fail{_RESET}]"
    if state == "building":
        return f"[{_YELLOW}building{_RESET}]"
    return "[pending]"


def _resultLine(galleryConf, index, fname, state, duration=None):
    """Return the permanent line of one example.

    The line is printed by the process that ran the example, so the terminal
    is detected from the gallery configuration here: inside the pipe capture
    of a worker ``sys.stdout`` is not a terminal anymore.
    """
    total = len(_examples[0]) if _examples else 0
    width = _nameWidth(galleryConf)
    name = _shorten(_displayName(index, fname), width)
    colour = _isLive(galleryConf) or sys.stdout.isatty()
    status = _statusText(state) if colour else f"[{state}]"
    line = f"[{index:3d}/{total}] {name:<{width}} ..."
    if duration is not None:
        return f"{line} {status} {duration:9.1f} s"
    return f"{line} {status}"


def _start(galleryConf, fname, index):
    """Report the start and redirect stdout and stderr into the example log.

    The messages of the C++ core go to stderr, so both descriptors have to be
    redirected. Otherwise a warning of one example appears on the console
    between the lines of other examples and cannot be attributed.
    """
    global _current

    _current = None
    _setStatus(galleryConf, index, "building")
    if not _isLive(galleryConf):
        print(_resultLine(galleryConf, index, fname, "building"), flush=True)

    savedFd = os.dup(1)
    savedErr = os.dup(2)
    try:
        logDir = _galleryDir(galleryConf) / "_gallery_logs"
        logDir.mkdir(parents=True, exist_ok=True)
        logFile = logDir / f"{Path(str(fname)).stem}.log"
        sink = os.open(logFile, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
        os.dup2(sink, 1)
        os.dup2(sink, 2)
        os.close(sink)
    except OSError:
        os.close(savedFd)
        os.close(savedErr)
        return

    _current = (savedFd, savedErr, fname, index, time.time())


def _finish(galleryConf):
    """Restore the descriptors and report the runtime of the example."""
    global _current

    if _current is None:
        return
    savedFd, savedErr, fname, index, started = _current
    _current = None
    duration = time.time() - started
    os.dup2(savedFd, 1)
    os.dup2(savedErr, 2)
    os.close(savedFd)
    os.close(savedErr)

    _setStatus(galleryConf, index, f"ok {duration:9.1f} s")
    _storeTime(galleryConf, fname, duration)
    if not _isLive(galleryConf):
        # While the table is drawn, the process that draws it reports the
        # finished examples, the output of a worker may be swallowed.
        print(
            _resultLine(galleryConf, index, fname, "ok", duration), flush=True
        )


def markFailed(galleryConf, fname):
    """Mark an example as failed, called when it raised an exception.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration.
    fname : str
        Name of the example that failed.
    """
    global _current

    if not _enabled or fname is None:
        return
    if _current is not None:
        # The example is aborted, so the "after" hook will not restore this.
        # Otherwise the abort message would end up in the example log.
        savedFd, savedErr = _current[0], _current[1]
        _current = None
        os.dup2(savedFd, 1)
        os.dup2(savedErr, 2)
        os.close(savedFd)
        os.close(savedErr)
    index = _slot(galleryConf, fname)
    if _readStatus(galleryConf).get(index, "").startswith("fail"):
        return
    _setStatus(galleryConf, index, "fail")
    if not _isLive(galleryConf):
        print(_resultLine(galleryConf, index, fname, "fail"), flush=True)


def reportProgress(galleryConf, fname, when="before"):
    """Report the progress of a gallery example.

    Used as a sphinx-gallery ``reset_modules`` hook with
    ``reset_modules_order = "both"``.

    Parameters
    ----------
    galleryConf : dict
        Sphinx-gallery configuration.
    fname : str
        Name of the example that is about to run or has been executed.
    when : str ["before"]
        Whether the hook is called before or after the example.
    """
    if not _enabled or fname is None or _isWarmed(galleryConf):
        return

    index = _slot(galleryConf, fname)
    if when == "before":
        _start(galleryConf, fname, index)
    else:
        _finish(galleryConf)


def _quietStatusIterator(app):
    """Drop the progress lines sphinx-gallery would print per example."""
    genRst = sys.modules.get("sphinx_gallery.gen_rst")
    if not _enabled or genRst is None:
        return
    if getattr(genRst.status_iterator, "__gimliQuiet__", False):
        return

    def statusIterator(iterable, *args, **kwargs):
        return iter(iterable)

    statusIterator.__gimliQuiet__ = True
    genRst.status_iterator = statusIterator


def _condenseWarnings():
    """Collect the warnings of the build, they belong to the gallery.

    The pool warnings of joblib come with the source path and two extra lines,
    which ended up between the example lines and could not be attributed. They
    are reported in one block when the examples are done.

    Returns
    -------
    The replaced hook, so the caller can restore it.
    """
    original = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None):
        _warnings.append(
            (category.__name__, str(message), Path(str(filename)).name, lineno)
        )

    warnings.showwarning = showwarning
    return original


def _startTable(app):
    """Start the printer that keeps the running examples at the end."""
    global _printer

    galleryConf = getattr(app.config, "sphinx_gallery_conf", None) or {}
    _quietStatusIterator(app)
    if not _enabled or not galleryConf:
        return
    app.gimliWarningHook = _condenseWarnings()
    # Markers of an aborted build must not decide about this one.
    _liveMarker(galleryConf).unlink(missing_ok=True)
    resetProgress(galleryConf)
    if not _live or not sys.stdout.isatty():
        return

    names = _orderedExamples(galleryConf)
    if not names:
        return

    printer = _TablePrinter(galleryConf, names)
    print(
        f"gallery: {len(names)} examples, finished ones as lines, "
        f"running ones as table, terminal "
        f"{printer.nameWidth + 36} of {printer.columns} columns",
        flush=True,
    )

    _statusDir(galleryConf).mkdir(parents=True, exist_ok=True)
    _liveMarker(galleryConf).write_text("live", encoding="utf-8")
    app.gimliPrinter = printer
    _printer = printer
    printer.capture()
    printer.start()


def _stopTable(app, exception):
    """Stop the printer and report where the example output is."""
    global _printer

    printer = getattr(app, "gimliPrinter", None)
    galleryConf = getattr(app.config, "sphinx_gallery_conf", None) or {}
    hook = getattr(app, "gimliWarningHook", None)
    if hook is not None:
        warnings.showwarning = hook
    if printer is not None:
        printer.release()
        _printer = None
        _liveMarker(galleryConf).unlink(missing_ok=True)
    _reportWarnings()

    logDir = _galleryDir(galleryConf) / "_gallery_logs"
    if not _enabled or not logDir.is_dir():
        return
    logs = sorted(logDir.glob("*.log"))
    if len(logs):
        print(f"\nConsole output of {len(logs)} examples in {logDir}")


def setup(app):
    """Register the extension.

    Parameters
    ----------
    app : sphinx.application.Sphinx
        The Sphinx application.

    Returns
    -------
    dict
        Extension metadata.
    """
    app.connect("builder-inited", _startTable, priority=100)
    app.connect("build-finished", _stopTable)
    return {
        "version": "1.0",
        "parallel_read_safe": False,
        "parallel_write_safe": True,
    }
