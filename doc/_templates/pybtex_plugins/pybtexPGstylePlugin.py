#!/usr/bin/env python3
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
#from pybtex.style.template import toplevel # ... and anything else needed
from pybtex.style.labels import BaseLabelStyle

class PGStyle(UnsrtStyle):
    name = 'pgstyle'
    default_name_style = 'lastfirst' # 'lastfirst' or 'plain'
    default_label_style = 'alpha' # 'number' or 'alpha'
    default_sorting_style = 'author_year_title' # 'none' or 'author_year_title'

    #def format_XXX(self, e):
        #template = toplevel [
            ## etc.
        #]
        #return template.format_data(e)

class Alpha(BaseLabelStyle):
    name = 'alpha'

    def format(self, entry):
        print('########  templates.pybtex.Alpha: #############')
        print(entry.__dict__)
        return str(entry.key)
