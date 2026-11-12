"""Create the publication list for the documentation from a bib file."""

import json
from pathlib import Path


def _entryToDict(entry):
    """Convert a bibtexparser 2 entry into a plain bibtex dictionary."""
    data = {"ID": entry.key, "ENTRYTYPE": entry.entry_type}
    for name, field in entry.fields_dict.items():
        data[name] = field.value
    return data


def parseBib(fname):
    """Read a bibtex file and return the entries sorted by year."""
    import bibtexparser

    if hasattr(bibtexparser, "parse_file"):  # bibtexparser 2.x
        library = bibtexparser.parse_file(fname)
        references = [_entryToDict(entry) for entry in library.entries]
    else:  # bibtexparser 1.x
        from bibtexparser.bparser import BibTexParser
        from bibtexparser.customization import convert_to_unicode

        parser = BibTexParser()
        parser.customization = convert_to_unicode
        with Path(fname).open() as bibfile:
            library = bibtexparser.load(bibfile, parser=parser)
        references = library.get_entry_list()

    references.sort(key=lambda x: x["year"], reverse=True)

    return references


def write_html(bibfile="gimliuses.bib"):
    """Return the publication list as JSON ready for the HTML template."""
    db = parseBib(bibfile)
    for entry in db:
        entry["author"] = entry["author"].replace(" and ", ", ")
        entry["author"] = entry["author"].replace("~", " ")
        if "journal" not in entry:
            entry["journal"] = entry.pop("booktitle", "JOURNAL_MISSING")
        entry["journal"] = f"<i>{entry['journal']}</i>"
        if "doi" in entry:
            link = f"https://doi.org/{entry['doi']}"
            string = (
                f"<a target='_blank' href={link} data-toggle='tooltip' "
                f"title='Go to {link}'><i class='ai ai-doi ai-2x'></i></a>"
            )
        else:
            string = ""
        entry["doi"] = string

    return json.dumps(db, sort_keys=True, indent=4)
