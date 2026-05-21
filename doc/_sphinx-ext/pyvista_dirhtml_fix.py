# Needed until https://github.com/pyvista/pyvista/issues/8688 is fixed
"""Replacement for pyvista.ext.viewer_directive that works under -b dirhtml.

The upstream directive builds the iframe `src` from the .rst source-tree depth
and reuses it as the URL from the rendered HTML page to `_static/viewer.html`.
That assumption only holds for the default `html` builder; `dirhtml` writes
each page as `<page>/index.html`, which is one level deeper, so the iframe URL
ends up too shallow and the viewer 404s. We override the directive and let
Sphinx's own URI machinery compute the path.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.util import logging
from sphinx.util.osutil import relative_uri
from trame_vtk.tools.vtksz2html import HTML_VIEWER_PATH

logger = logging.getLogger(__name__)


class OfflineViewerDirective(Directive):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    has_content = True

    def run(self):
        env = self.state.document.settings.env
        builder = env.app.builder
        source_dir = Path(env.app.srcdir)
        output_dir = Path(env.app.outdir)
        build_dir = output_dir.parent

        source_file = (
            Path(self.state.document.current_source).parent / self.arguments[0]
        ).absolute().resolve()
        if not source_file.is_file():
            logger.warning(f"Source file {source_file} does not exist.")
            return []

        # Copy the viewer HTML into _static once.
        static_path = output_dir / "_static"
        static_path.mkdir(exist_ok=True)
        viewer_name = Path(HTML_VIEWER_PATH).name
        if not (static_path / viewer_name).exists():
            shutil.copy(HTML_VIEWER_PATH, static_path)

        # Mirror the source layout under _images/ to avoid collisions between
        # vtksz files generated from different pages.
        if source_file.is_relative_to(build_dir):
            dest_partial_path = source_file.parent.relative_to(build_dir)
        elif source_file.is_relative_to(source_dir):
            dest_partial_path = source_file.parent.relative_to(source_dir)
        else:
            logger.warning(
                f"Source file {source_file} is not inside the build or source "
                "directory; cannot place asset."
            )
            return []

        dest_path = output_dir / "_images" / dest_partial_path
        dest_path.mkdir(parents=True, exist_ok=True)
        dest_file = (dest_path / source_file.name).resolve()
        if source_file != dest_file:
            try:
                shutil.copy(source_file, dest_file)
            except Exception as e:  # noqa: BLE001
                logger.warning(f"Failed to copy {source_file} to {dest_file}: {e}")

        # Builder-aware URIs. `get_target_uri(docname)` gives the right output
        # path for the current page in any builder (e.g. `foo/bar.html` for
        # html, `foo/bar/` for dirhtml), so `relative_uri` produces a link that
        # works in both.
        page_uri = builder.get_target_uri(env.docname)
        viewer_uri_in_output = f"_static/{viewer_name}"
        asset_uri_in_output = (
            Path("_images") / dest_partial_path / source_file.name
        ).as_posix()

        rel_viewer_path = relative_uri(page_uri, viewer_uri_in_output)
        # The viewer resolves `fileURL` relative to itself, so this path is
        # builder-independent.
        rel_asset_path = relative_uri(viewer_uri_in_output, asset_uri_in_output)

        html = (
            f"<iframe src='{rel_viewer_path}?fileURL={rel_asset_path}' "
            "width='100%' height='400px' frameborder='0'></iframe>"
        )
        return [nodes.raw("", html, format="html")]


def setup(app):
    # Override pyvista.ext.viewer_directive's registration. Sphinx's
    # add_directive replaces any previous directive of the same name, so as
    # long as this extension is listed after `pyvista.ext.viewer_directive`
    # in conf.py, our class wins.
    app.add_directive("offlineviewer", OfflineViewerDirective, override=True)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
