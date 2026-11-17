"""Execute all gallery examples in a single parallel pool.

sphinx-gallery executes the examples of a gallery section in its own
``joblib.Parallel`` block and waits for that block before the next section is
started, because the section index needs the results. With many small sections
that leaves most of the workers idle at every section boundary.

Here all examples of all sections are executed once up front in one pool. The
regular gallery run then finds every example up to date (same md5) and only
assembles the ``index.rst`` files, which takes a few seconds per section.

Set ``GIMLI_DOC_GLOBAL_PARALLEL=0`` to get the sphinx-gallery behaviour back.
"""

import os
import sys
from pathlib import Path

__all__ = ["setup"]

_enabled = os.environ.get("GIMLI_DOC_GLOBAL_PARALLEL", "1") != "0"


def _folders(galleryConf, srcRoot, getSubsections, prepareDirs):
    """Return the (source, target) folder pairs of all gallery sections."""
    folders = []
    for exampleDir, galleryDir in prepareDirs(galleryConf, srcRoot):
        examplePath = srcRoot / str(exampleDir)
        galleryPath = srcRoot / str(galleryDir)
        folders.append((examplePath, galleryPath))
        for subsection in getSubsections(srcRoot, examplePath, galleryConf):
            name = Path(str(subsection)).name
            folders.append((examplePath / name, galleryPath / name))
    return folders


def _collect(galleryConf, folders, collectFiles, getCallables, splitParallel):
    """Return the examples to execute as (name, source dir, target dir).

    Examples that opt out of parallel execution are returned separately, they
    have to run alone in this process like sphinx-gallery does it.
    """
    serial = []
    parallel = []
    for srcDir, targetDir in folders:
        (targetDir / "images" / "thumb").mkdir(parents=True, exist_ok=True)
        listdir = collectFiles([srcDir], galleryConf)
        (sortKey,) = getCallables(
            galleryConf, "within_subsection_order", srcDir
        )
        sortedListdir = sorted(listdir, key=sortKey)
        serialList, parallelList = splitParallel(
            sortedListdir, srcDir, galleryConf
        )
        serial.extend((name, srcDir, targetDir) for name in serialList)
        parallel.extend((name, srcDir, targetDir) for name in parallelList)
    return serial, parallel


def _markWarmed(galleryConf):
    """Tell the progress table that the examples have been executed."""
    try:
        from galleryprogress import markWarmed
    except ImportError:
        return
    markWarmed(galleryConf)


def _flushProgress():
    """Report the examples that finished since the last refresh.

    The table is refreshed from a polling thread, so without this the line of
    the last example would appear after the output that sphinx writes next.
    """
    try:
        from galleryprogress import flushProgress
    except ImportError:
        return
    flushProgress()


def _finishProgress():
    """Remove the status block, the examples have been executed."""
    try:
        from galleryprogress import finishProgress
    except ImportError:
        return
    finishProgress()


def _resetProgress(galleryConf):
    """Forget the states of the previous run."""
    try:
        from galleryprogress import resetProgress
    except ImportError:
        return
    resetProgress(galleryConf)


def _exampleTimes(galleryConf):
    """Return the cached runtime per example name."""
    try:
        from galleryprogress import exampleTimes
    except ImportError:
        return {}
    return exampleTimes(galleryConf)


def _slowestFirst(galleryConf, tasks):
    """Sort the examples by runtime, slowest first.

    The scheduling of the last run is used, so that the examples that need
    most of the time are started right away and do not end up alone at the end.
    Unknown examples, e.g. new ones, keep their position at the end.
    """
    times = _exampleTimes(galleryConf)
    tasks.sort(key=lambda task: -times.get(Path(str(task[0])).name, 0.0))
    return tasks


def _warmup(app):
    """Execute every gallery example once, in one pool of workers."""
    galleryConf = getattr(app.config, "sphinx_gallery_conf", None) or {}
    if not _enabled or not galleryConf or not galleryConf.get("parallel"):
        return

    try:
        from joblib import Parallel, delayed
        from sphinx_gallery.gen_gallery import (
            _prepare_sphx_glr_dirs,
            get_subsections,
        )
        from sphinx_gallery.gen_rst import (
            _get_callables,
            _split_parallel,
            generate_file_rst,
        )
        from sphinx_gallery.utils import _collect_gallery_files
    except ImportError as err:
        print(f"global gallery parallelisation disabled: {err}")
        return

    srcRoot = Path(app.builder.srcdir)
    try:
        folders = _folders(
            galleryConf, srcRoot, get_subsections, _prepare_sphx_glr_dirs
        )
        serial, tasks = _collect(
            galleryConf,
            folders,
            _collect_gallery_files,
            _get_callables,
            _split_parallel,
        )
    except (AttributeError, OSError, TypeError) as err:
        print(f"global gallery parallelisation disabled: {err}")
        return

    if not serial and not tasks:
        return

    _slowestFirst(galleryConf, tasks)
    _resetProgress(galleryConf)

    if not sys.stdout.isatty():
        print(
            f"gallery: {len(serial) + len(tasks)} examples in "
            f"{len(folders)} sections executed in one pool of "
            f"{galleryConf['parallel']}",
            flush=True,
        )

    for name, srcDir, targetDir in serial:
        generate_file_rst(name, targetDir, srcDir, galleryConf)

    if tasks:
        pool = Parallel(
            n_jobs=galleryConf["parallel"],
            pre_dispatch="2 * n_jobs",
            batch_size=1,
            backend="loky",
        )
        pool(
            delayed(generate_file_rst)(name, targetDir, srcDir, galleryConf)
            for name, srcDir, targetDir in tasks
        )

    _flushProgress()
    _finishProgress()
    _markWarmed(galleryConf)


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
    # After the progress table is started, but before sphinx-gallery runs.
    app.connect("builder-inited", _warmup, priority=200)
    return {
        "version": "1.0",
        "parallel_read_safe": False,
        "parallel_write_safe": True,
    }
