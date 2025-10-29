#!/usr/bin/env python3
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
#from pybtex.style.template import toplevel # ... and anything else needed
from pybtex.style.labels.alpha import LabelStyle as AlphaLabelStyle


class KeyLabelStyle(AlphaLabelStyle):
    def format_label(self, entry):
        label = entry.key
        return str(label)


class PGStyle(UnsrtStyle):
    name = 'pgstyle'
    default_sorting_style = 'author_year_title'
    default_name_style = 'lastfirst' # 'lastfirst' or 'plain'
    default_label_style = 'alpha' # 'number' or 'alpha'

    def __init__(self, *args, **kwargs):
        super(PGStyle, self).__init__(*args, **kwargs)
        self.label_style = KeyLabelStyle()
        self.format_labels = self.label_style.format_labels