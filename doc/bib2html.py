import json
import re

import bibtexparser

# bibtexparser 2.x dropped the 1.x API (load/BibTexParser/customization),
# so support both as long as either version may be installed.
BIBTEXPARSER_V2 = not hasattr(bibtexparser, "load")


def _parse_bib_v1(fname):
    """ Read bibtex file with the bibtexparser 1.x API. """
    from bibtexparser.bparser import BibTexParser
    from bibtexparser.customization import convert_to_unicode

    with open(fname) as bibfile:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bp = bibtexparser.load(bibfile, parser=parser)

    return bp.get_entry_list()


def _latex_decoder_v2():
    """ Build a LaTeX decoder for bibtexparser 2.x.

    Same as the one LatexDecodingMiddleware creates by default, but keeps
    stray ampersands intact. Some entries contain a plain (or even HTML
    escaped) "&" instead of "\\&", which pylatexenc would otherwise swallow
    as a tabular alignment character.
    """
    from pylatexenc.latex2text import (LatexNodes2Text, MacroTextSpec,
                                       SpecialsTextSpec,
                                       get_default_latex_context_db)

    context = get_default_latex_context_db()
    context.add_context_category(
        "pygimli-bib-context",
        prepend=True,
        macros=[MacroTextSpec("url", simplify_repl="%s")],  # no '< ... >'
        specials=[SpecialsTextSpec("&", "&")],
    )

    return LatexNodes2Text(latex_context=context, math_mode="verbatim")


def _parse_bib_v2(fname):
    """ Read bibtex file with the bibtexparser 2.x API. """
    from bibtexparser.middlewares import LatexDecodingMiddleware

    layers = [LatexDecodingMiddleware(decoder=_latex_decoder_v2())]
    library = bibtexparser.parse_file(fname, append_middleware=layers)

    entries = []
    for entry in library.entries:
        fields = {field.key: field.value for field in entry.fields}
        fields["ID"] = entry.key
        fields["ENTRYTYPE"] = entry.entry_type
        entries.append(fields)

    return entries


def parse_bib(fname):
    """ Read bibtex file and sort by year. """
    if BIBTEXPARSER_V2:
        references = _parse_bib_v2(fname)
    else:
        references = _parse_bib_v1(fname)

    references.sort(key=lambda x: x["year"], reverse=True)

    return references


def write_html(bibfile="gimliuses.bib"):
    db = parse_bib(bibfile)
    for entry in db:
        # "~" may already have been decoded into a non-breaking space
        entry["author"] = entry["author"].replace("~", " ").replace("\u00a0", " ")
        entry["author"] = re.sub(r"\s+and\s+", ", ", entry["author"])
        if not "journal" in entry:
            entry["journal"] = entry.pop("booktitle", "JOURNAL_MISSING")
        entry["journal"] = "<i>%s</i>" % entry["journal"]
        if not "doi" in entry:
            string = ""
        else:
            doi = entry["doi"]
            link = "https://doi.org/%s" % doi
            string = (
                "<a target='_blank' href=%s data-toggle='tooltip' title='Go to %s'><i class='ai ai-doi ai-2x'></i></a>"
                % (link, link)
            )
        entry["doi"] = string

    return json.dumps(db, sort_keys=True, indent=4)
