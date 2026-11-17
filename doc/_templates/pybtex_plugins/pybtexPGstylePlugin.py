#!/usr/bin/env python3
"""Pybtex style plugin for pyGIMLi documentation builds.

Registers a citation style that labels entries with their bibtex key
and sorts them by author, year and title.
"""

from pybtex.style.formatting.unsrt import Style as UnsrtStyle

# from pybtex.style.template import toplevel # ... and anything else needed
from pybtex.style.labels.alpha import LabelStyle as AlphaLabelStyle


class KeyLabelStyle(AlphaLabelStyle):
    """Label style that uses the bibtex key of an entry as its label."""

    def format_label(self, entry):
        """Return the bibtex key of the given entry as label."""
        label = entry.key
        return str(label)


class PGStyle(UnsrtStyle):
    """Unsrt style with key-based labels for pyGIMLi citations."""

    name = "pgstyle"
    default_sorting_style = "author_year_title"
    default_name_style = "lastfirst"  # 'lastfirst' or 'plain'
    default_label_style = "alpha"  # 'number' or 'alpha'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.label_style = KeyLabelStyle()
        self.format_labels = self.label_style.format_labels
