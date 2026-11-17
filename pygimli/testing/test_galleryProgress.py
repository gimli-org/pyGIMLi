"""Regression tests for the documentation gallery progress output.

The extension lives in ``doc/_sphinx-ext`` and is not part of the installed
package, so the tests are skipped when the source tree is missing.

Covered failures:

* a warning of the C++ core (gmsh) is written to stderr, so redirecting only
  stdout let it appear on the console between the lines of other examples;
* the rows of the block were separated by a line feed only, which does not
  reset the cursor column and looked like duplicated blocks in the terminal;
* a finished example has to leave the block as a permanent line, exactly once.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

EXT_DIR = Path(__file__).resolve().parents[2] / "doc" / "_sphinx-ext"
if EXT_DIR.is_dir():
    sys.path.insert(0, str(EXT_DIR))

try:
    import galleryprogress
except ImportError:  # pragma: no cover - installed package without the docs
    galleryprogress = None


@unittest.skipIf(
    galleryprogress is None, "documentation extensions not available"
)
class TestGalleryProgress(unittest.TestCase):
    """Progress output of the sphinx gallery."""

    def setUp(self):
        """Create a fake gallery and remember the module state."""
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.examples = galleryprogress._examples
        self.current = galleryprogress._current
        galleryprogress._examples = None
        galleryprogress._current = None

        self.exampleDir = self.root / "examples"
        self.exampleDir.mkdir()
        for index in range(1, 4):
            (self.exampleDir / f"plot_{index:02d}_demo.py").write_text(
                "pass\n", encoding="utf-8"
            )
        (self.root / "auto_examples" / "_gallery_status").mkdir(parents=True)
        self.conf = {
            "examples_dirs": ["examples"],
            "gallery_dirs": ["auto_examples"],
            "filename_pattern": "/plot_",
            "src_dir": str(self.root),
        }

    def tearDown(self):
        """Restore the module state and remove the temporary files."""
        galleryprogress._examples = self.examples
        galleryprogress._current = self.current
        self.directory.cleanup()

    def printer(self):
        """Return a printer that writes into a file instead of a terminal."""
        printer = galleryprogress._TablePrinter(
            self.conf, galleryprogress._orderedExamples(self.conf)
        )
        fd = os.open(
            self.root / "screen.txt", os.O_RDWR | os.O_CREAT | os.O_TRUNC
        )
        self.addCleanup(os.close, fd)
        printer.terminalFd = fd
        printer.columns = 100
        printer.room = 20
        printer.nameWidth = 60
        return printer, fd

    @staticmethod
    def content(fd):
        """Return everything that was written to the printer file."""
        os.lseek(fd, 0, os.SEEK_SET)
        return os.read(fd, 1 << 20).decode("utf-8", "replace")

    def test_stderr_of_an_example_is_kept_in_its_log(self):
        """A warning on stderr belongs to the example that caused it."""
        galleryprogress._start(self.conf, "examples/plot_01_demo.py", 1)
        try:
            os.write(1, b"stdout of the example\n")
            os.write(2, b"Warning : Unknown entity of dimension 2\n")
        finally:
            galleryprogress._finish(self.conf)

        log = self.root / "auto_examples" / "_gallery_logs" / "plot_01_demo.log"
        text = log.read_text(encoding="utf-8")
        self.assertIn("stdout of the example", text)
        self.assertIn("Warning : Unknown entity", text)

    def test_failed_example_restores_the_descriptors(self):
        """An aborted example must not keep stdout and stderr redirected."""
        galleryprogress._start(self.conf, "examples/plot_02_demo.py", 2)
        galleryprogress.markFailed(self.conf, "examples/plot_02_demo.py")

        self.assertIsNone(galleryprogress._current)
        status = self.root / "auto_examples" / "_gallery_status" / "2"
        self.assertTrue(status.read_text(encoding="utf-8").startswith("fail"))

    def test_table_rows_are_written_with_crlf(self):
        r"""Rows need a carriage return, a line feed keeps the column.

        Without the explicit ``\r`` the rows of the block start where the
        previous row ended, which looked like duplicated blocks.
        """
        printer, fd = self.printer()
        for index in (1, 2):
            galleryprogress._setStatus(self.conf, index, "building")
        printer.render()

        text = self.content(fd)
        self.assertIn("[  1/3] examples/plot_01_demo.py", text)
        self.assertIn("[  2/3] examples/plot_02_demo.py", text)
        self.assertIn("\r\n", text)
        self.assertEqual(text.count("\x1b[33m"), 2)
        self.assertIn("0/3 done, 2 running", text)

    def test_finished_example_is_written_once(self):
        """A finished example leaves the block as one permanent line."""
        printer, fd = self.printer()
        galleryprogress._setStatus(self.conf, 1, "building")
        printer.render()
        galleryprogress._setStatus(self.conf, 1, "ok      12.3 s")
        printer.render()
        printer.render()

        text = self.content(fd)
        self.assertEqual(text.count("\x1b[32mok"), 1)
        self.assertIn("12.3 s", text)
        tail = text[text.index("\x1b[32mok") :]
        self.assertNotIn("building", tail)
        self.assertIn("1/3 done, 0 running", tail)

    def test_erase_moves_up_over_the_whole_block(self):
        """The erase asks for one row less than the block has."""
        printer, fd = self.printer()
        for index in (1, 2, 3):
            galleryprogress._setStatus(self.conf, index, "building")
        printer.render()
        printer.erase()

        text = self.content(fd)
        # Three running rows plus the footer: up three rows, clear below.
        self.assertIn("\x1b[3A\r\x1b[J", text)

    def test_reset_forgets_the_warmed_marker(self):
        """A marker of the previous run must not switch the hooks off.

        ``_warmed`` survives in the gallery directory, and a build that finds
        it skips the redirects of the examples, so their output floods the
        terminal again.
        """
        galleryprogress._setStatus(self.conf, 1, "ok      1.0 s")
        galleryprogress.markWarmed(self.conf)
        galleryprogress.resetProgress(self.conf)

        self.assertFalse(galleryprogress._isWarmed(self.conf))
        self.assertEqual(galleryprogress._readStatus(self.conf), {})

    def test_reset_keeps_the_live_marker(self):
        """The warm-up resets the states while the table is already running."""
        galleryprogress._liveMarker(self.conf).write_text(
            "live", encoding="utf-8"
        )
        galleryprogress.resetProgress(self.conf)

        self.assertTrue(galleryprogress._isLive(self.conf))

    def test_hooks_redirect_the_output_after_a_reset(self):
        """After the reset the next example writes into its log again."""
        galleryprogress.markWarmed(self.conf)
        galleryprogress.resetProgress(self.conf)

        galleryprogress.reportProgress(
            self.conf, "examples/plot_03_demo.py", "before"
        )
        try:
            self.assertIsNotNone(galleryprogress._current)
            os.write(1, b"output of the third example\n")
        finally:
            galleryprogress._finish(self.conf)

        log = self.root / "auto_examples" / "_gallery_logs" / "plot_03_demo.log"
        self.assertIn(
            "output of the third example", log.read_text(encoding="utf-8")
        )


if __name__ == "__main__":
    unittest.main()
