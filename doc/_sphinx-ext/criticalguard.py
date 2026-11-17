"""Abort the documentation build when a CRITICAL message is logged.

The C++ core reports unrecoverable problems with ``log(Critical, ...)``
(for example a forward response requested for a model with non-positive
resistivity). In contrast to the Python function :func:`pygimli.critical`,
which raises an exception, the C++ variant only writes a log message and
returns. An example with such a flaw is therefore still rendered and the
documentation build passes silently.

Sphinx-gallery executes the examples inside the Sphinx process, so the
messages are recorded here and turned into an exception after each executed
code block. Together with ``abort_on_example_error`` this fails the example
and thus the whole build with a readable message, and it keeps the faulty
example from being cached as successfully executed.

CRITICAL messages that are logged outside of an example (e.g. while the
configuration or the API documentation is imported) abort the build as well.

With a parallel gallery (``sphinx_gallery_conf["parallel"]``) every example is
executed in its own worker process, so the recorder is installed lazily in
each of them (see :func:`_ensureRecorder`). Messages of the build itself are
then reported when the build finishes.
"""

import contextlib
import logging
import os
import sys
import traceback
from pathlib import Path

__all__ = ["abortOnCriticalLog", "criticalLogScraper", "startExample"]

_loggerNames = ("core", "pygimli")
_criticalMessages = []
_buildMessages = []


class _CriticalRecorder(logging.Handler):
    """Collect CRITICAL messages of the pyGIMLi loggers."""

    def emit(self, record):
        """Store *record* if it is a pyGIMLi CRITICAL message."""
        name = record.name.split(".")[0].lower()
        if record.levelno >= logging.CRITICAL and name in _loggerNames:
            _criticalMessages.append(
                f"{record.name}: {record.getMessage()} "
                f"(logged at {record.pathname}:{record.lineno})"
            )


def _ensureRecorder():
    """Install the recorder unless this process has one already.

    With a parallel gallery (``sphinx_gallery_conf["parallel"]``) every
    example runs in a worker process, which only imports this module when the
    scraper is unpickled. Calling this from the hooks covers the case that the
    module was imported without running the module level installation.
    """
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, _CriticalRecorder):
            return
    root.addHandler(_CriticalRecorder())


def _abort(messages, context):
    """Raise a RuntimeError listing the recorded CRITICAL messages."""
    records = "\n".join(f"  {msg}" for msg in messages)
    messages.clear()
    raise RuntimeError(
        f"CRITICAL message(s) logged by pyGIMLi {context}, aborting the "
        f"documentation build:\n{records}"
    )


def _relativeSource(srcFile, galleryConf):
    """Return *srcFile* relative to the gallery source directory."""
    srcDir = str(galleryConf.get("src_dir", ""))
    srcFile = str(srcFile)
    if srcDir and srcFile.startswith(srcDir):
        return os.path.relpath(srcFile, srcDir)
    return srcFile


def _exampleLocation(exc, srcFile):
    """Return the innermost traceback position inside the example.

    Returns
    -------
    str
        Empty string if the example itself does not show up in the traceback.
    """
    base = Path(str(srcFile)).name
    for frame in reversed(traceback.extract_tb(exc.__traceback__)):
        if Path(frame.filename).name != base:
            continue
        location = f"line {frame.lineno}"
        if frame.line:
            location += f": {frame.line.strip()}"
        return location
    return ""


def _ensureExampleNames():
    """Report the example that raised an exception.

    Sphinx-gallery catches the exception of a code block and logs which file
    failed, but that log output goes through its own stdout/stderr tee and is
    therefore lost in the build log. Only the bare exception reaches Sphinx and
    aborts the build without telling which example it came from. This wraps the
    block execution so the message names the example and the line instead.
    """
    genRst = sys.modules.get("sphinx_gallery.gen_rst")
    if genRst is None:
        return

    original = genRst.execute_code_block
    if getattr(original, "__gimliWrapped__", False):
        return

    def executeCodeBlock(
        compiler, block, exampleGlobals, scriptVars, galleryConf, fileConf
    ):
        try:
            return original(
                compiler,
                block,
                exampleGlobals,
                scriptVars,
                galleryConf,
                fileConf,
            )
        except BaseException as exc:
            srcFile = scriptVars.get("src_file", "<unknown example>")
            with contextlib.suppress(ImportError, AttributeError, OSError):
                # Show the example as failed in the progress table.
                import galleryprogress

                galleryprogress.markFailed(galleryConf, srcFile)
            where = _relativeSource(srcFile, galleryConf)
            lineno = getattr(block, "lineno", None)
            if lineno:
                where = f"{where} (code block starting at line {lineno})"
            location = _exampleLocation(exc, srcFile)
            if location:
                where = f"{where}, {location}"
            raise RuntimeError(
                f"{exc}\n\nRaised while running {where}, aborting the "
                f"documentation build."
            ) from exc

    executeCodeBlock.__gimliWrapped__ = True
    genRst.execute_code_block = executeCodeBlock


def startExample(gallery_conf, fname):
    """Split off CRITICAL messages that are not part of an example.

    Used as a sphinx-gallery ``reset_modules`` hook, which is called before
    every executed example, in the Sphinx process as well as in every gallery
    worker process. Messages logged by the build itself (e.g. while importing
    the API documentation) must not be attributed to whichever example happens
    to run next, so they are kept separately.

    Parameters
    ----------
    gallery_conf : dict
        Sphinx-gallery configuration, unused.
    fname : str
        Name of the example that is about to run, unused.
    """
    _ensureRecorder()
    _ensureExampleNames()
    _buildMessages.extend(_criticalMessages)
    _criticalMessages.clear()


def criticalLogScraper(block, block_vars, gallery_conf):
    """Abort the build if a CRITICAL message was logged.

    Used as a sphinx-gallery image scraper, which runs after every executed
    code block.

    Parameters
    ----------
    block : sphinx_gallery.py_source_parser.Block
        Executed code block, only used for its line number.
    block_vars : dict
        Sphinx-gallery runtime variables of the current example, used for the
        source file name.
    gallery_conf : dict
        Sphinx-gallery configuration, used for the source directory.

    Returns
    -------
    str
        Always empty, no images are produced.

    Raises
    ------
    RuntimeError
        If a CRITICAL message was logged since the last check.
    """
    _ensureRecorder()

    if _buildMessages:
        _abort(_buildMessages, "outside of the gallery examples")

    if not _criticalMessages:
        return ""

    srcFile = _relativeSource(
        block_vars.get("src_file", "<unknown example>"), gallery_conf
    )

    lineno = getattr(block, "lineno", None)
    if lineno:
        where = f"{srcFile} (code block starting at line {lineno})"
    else:
        where = srcFile

    _abort(_criticalMessages, f"while running {where}")
    return ""


def abortOnCriticalLog(app, exception):
    """Fail the build if CRITICAL messages were logged without an example.

    Connected to ``build-finished``, which is the only chance to react when no
    gallery example was executed at all (e.g. with ``plot_gallery=0``).

    Parameters
    ----------
    app : sphinx.application.Sphinx
        Sphinx application, unused.
    exception : BaseException or None
        Exception of a failed build, then it is already reported.
    """
    if exception is not None:
        return

    if _buildMessages:
        _abort(_buildMessages, "outside of the gallery examples")
    if _criticalMessages:
        _abort(_criticalMessages, "outside of the gallery examples")


# The recorder must be active for every example code block, so it is
# installed as soon as this module is imported by the documentation config.
_ensureRecorder()


def setup(app):
    """Register the extension.

    Parameters
    ----------
    app : sphinx.application.Sphinx
        Sphinx application, unused.

    Returns
    -------
    dict
        Extension metadata.
    """
    app.connect("build-finished", abortOnCriticalLog)
    return {
        "version": "1.0",
        "parallel_read_safe": False,
        "parallel_write_safe": True,
    }
