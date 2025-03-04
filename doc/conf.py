# -*- coding: utf-8 -*-
"""
pyGIMLi sphinx configuration file.
"""

import warnings

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message="Matplotlib is currently using agg, which is a"
    " non-GUI backend, so cannot show the figure.",
)

import random
import datetime
import os
import re
import sys
from os import path
from os.path import join

sys.path.insert(0, os.path.abspath("."))

import numpy as np

# for doc rendering on headless machines (jenkins server)
import matplotlib

matplotlib.use("Agg")
import pkg_resources
import sphinx

import pygimli as pg


from sidebar_gallery import make_gallery

try:
    # from _build.doc.conf_environment import *
    from conf_environment import *

    pg.boxprint("Building documentation out-of-source. Good.")
    in_source = False
    print("DOXY_BUILD_DIR", DOXY_BUILD_DIR)
except ImportError:
    TRUNK_PATH = ".."
    SPHINXDOC_PATH = "."
    DOC_BUILD_DIR = ""
    DOXY_BUILD_DIR = ""
    in_source = True
    pg.boxprint("Building documentation in-source. Don't forget to make clean.")

sys.path.append(os.path.abspath(SPHINXDOC_PATH))
sys.path.append(os.path.abspath(join(SPHINXDOC_PATH, "_sphinx-ext")))

# The following line is necessary for the Tools section
sys.path.append(os.path.abspath(join(TRUNK_PATH, "pygimli")))

# -- General configuration ----------------------------------------------------

# MPL configuration in API docs, tutorials and examples
plot_rcparams = {"savefig.bbox": "tight"}

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "1.8"  # due to napoleon

# Check for external sphinx extensions
deps = [
    "sphinxcontrib-programoutput",
    "sphinxcontrib-bibtex",
    "sphinxcontrib-doxylink"
]

# check for p.version too
modules = [p.project_name for p in pkg_resources.working_set]

pg._b(modules)
req = []
for dep in deps:
    if dep not in modules:
        req.append(dep)
if req:
    msg = (
        "Sorry, there are missing dependencies to build the docs.\n"
        + "Try: sudo pip install %s.\n" % (" ".join(req))
        + "Or install all dependencies with: pip install -r dev_requirements.txt\n"
        + "You can install them all in userspace by adding the --user flag."
    )
    print((pkg_resources.working_set))
    raise ImportError(msg)

# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named "sphinx.ext.*")
# or your custom ones.
extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.imgconverter",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
    "matplotlib.sphinxext.plot_directive",
    "srclinks",
    "sphinxcontrib.doxylink",
    "sphinx_design",
    "sphinx_tippy"
    # "sphinxcontrib.spelling"
    #'sphinx.ext.pngmath',   # for breath
    #'sphinx.ext.todo',      # for breath
    #'breathe',              # doxgen to sphinx api docu
]

extensions += [dep.replace("-", ".") for dep in deps]


# Sphinx-gallery settings
try:
    import sphinx_gallery
    from sphinx_gallery.sorting import FileNameSortKey

    extensions += ["sphinx_gallery.gen_gallery"]

    def reset_mpl(gallery_conf, fname):
        import matplotlib

        matplotlib.rcParams.update(plot_rcparams)

    # Setup automatic gallery generation
    sphinx_gallery_conf = {
        "examples_dirs": [
            join(SPHINXDOC_PATH, "examples"),
            join(SPHINXDOC_PATH, "tutorials"),
        ],
        "gallery_dirs": ["_examples_auto", "_tutorials_auto"],
        # "reference_url": {
        #     "pygimli": "https://pygimli.org",
        #     "python": "https://docs.python.org/dev",
        #     "numpy": "https://numpy.org/doc/stable",
        #     "scipy": "https://docs.scipy.org/doc/scipy/reference",
        #     "matplotlib": "https://matplotlib.org/stable",
        # },
        # Don"t report time of fast scripts (< 10 sec)
        "min_reported_time": 10,
        # path where to store your example linker templates
        "backreferences_dir": "pygimliapi" + os.path.sep + "_generated",
        # Your documented modules. You can use a string or a list of strings
        "doc_module": "pygimli",
        # Sort gallery example by file name instead of number of lines (default)
        "within_subsection_order": FileNameSortKey,
        "remove_config_comments": True,
        # Only parse filenames starting with plot_
        "filename_pattern": "/plot_",
        "first_notebook_cell": ("# Checkout www.pygimli.org for more examples"),
        "reset_modules": (reset_mpl),
        # Avoid representation of mpl axis, LineCollections, etc.
        "ignore_repr_types": r"matplotlib[text, axes, collections]",
    }

    pyvista = pg.optImport("pyvista", "build the gallery with 3D visualizations")
    if pyvista:
        from pyvista.plotting.utilities.sphinx_gallery import DynamicScraper

        # necessary when building the sphinx gallery
        os.environ["PYVISTA_BUILDING_GALLERY"] = "true"
        pyvista.BUILDING_GALLERY = True
        os.environ["PYVISTA_OFF_SCREEN"] = "true"
        pyvista.OFF_SCREEN = True
        pyvista.set_plot_theme("document")
        pyvista.global_theme.window_size = np.array([1024, 768]) * 2
        pyvista.global_theme.font.size = 22
        pyvista.global_theme.font.label_size = 22
        pyvista.global_theme.font.title_size = 22
        pyvista.global_theme.return_cpos = False
        extensions += ["pyvista.ext.viewer_directive"]
        sphinx_gallery_conf["image_scrapers"] = (DynamicScraper(), "matplotlib")

except ImportError:
    err = """
    The sphinx_gallery extension is not installed, so the tutorials
    and examples will not be rendered and additional warnings will occur
    due to missing references.

    Install sphinx_gallery via:

        sudo pip3 install sphinx-gallery
    """
    pg.warn(err)


intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

autosectionlabel_prefix_document = True
autoclass_content = "class"
autosummary_generate = True
autosummary_generate_overwrite = False
autosummary_imported_members = True
autodoc_mock_imports = [
    "os",
    "os.path" "sys",
    "locale",
    "numpy",
    "matplotlib",
    "matplotlib.pyplot",
    "pyvista",
    "pyqt5",
]

autodoc_default_options = {
    "imported-members": True,
    # "special-members": "__init__",
    "undoc-members": False,
    "show-inheritance": True,
}

# Add any paths that contain templates here, relative to this directory.
templates_path = [join(SPHINXDOC_PATH, "_templates"), join(DOC_BUILD_DIR, "_templates")]

# MPL plot directive settings
plot_formats = [("png", 96)]
if not in_source:
    plot_formats.append(("pdf", 96))
plot_include_source = True
plot_html_show_source_link = False
plot_apply_rcparams = True  # if context option is used

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
source_encoding = "utf-8-sig"

# The master toctree document.
master_doc = "user-guide/index"

# General information about the project.
project = "pyGIMLi"
year = datetime.date.today().year
copyright = f"{year} - pyGIMLi Development Team"

# The version info for the project you"re documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = pg.__version__

rst_epilog = """
.. |version| replace:: pyGIMLi {versionnum}
""".format(versionnum=version)

# The full version, including alpha/beta/rc tags.
release = pg.__version__

release = release.replace("_", "\\_")

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ""
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = "%B %d, %Y"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = [
    "_build",
    "_sphinx-ext",
    "_templates",
    "tmp",
    "examples",
    "tutorials",
    ".venv",
    "../.venv",
]

# The reST default role (used for this markup: `text`) to use for all documents
# default_role = None

# If true, "()" will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# -- Options for HTML output --------------------------------------------------

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = [join(SPHINXDOC_PATH, "_themes")]

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "pydata_sphinx_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "logo": {
        "text": "py<b>GIMLi</b>",
    },
    "secondary_sidebar_items": ["page-toc", "improve-this-page"],
    "footer_start": ["footer_start"],
    "footer_end": ["footer_end"],
    "pygment_light_style": "friendly",
    "pygment_dark_style": "native",
}

# Temp: SEG announcement
import datetime
today = datetime.datetime.now()
webinar = datetime.datetime(2024, 3, 19)
if today < webinar:
    html_theme_options["announcement"] = "There will be a webinar on pyGIMLi hosted by SEG on March 19, 2024 at 4 pm CET. Register for free <a href='https://seg.org/calendar_events/open-source-software-webinar-pygimli/', target='_blank'>here</a>."

html_css_files = [
    "css/custom.css",
]

html_js_files = [
    "js/jquery-3.7.1.min.js",
]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = "pyGIMLi - Geophysical Inversion and Modelling Library"

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "pyGIMLi"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = join(SPHINXDOC_PATH, "_static/gimli_logo_simple.svg")

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = join(SPHINXDOC_PATH, "_static/favicon.ico")

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [join(SPHINXDOC_PATH, "_static")]

# If not "", a "Last updated on:" timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%b %d, %Y"  # + " with " + version

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    "index": [],
}

# Additional templates that should be rendered to pages, maps page names to
# template names.
html_additional_pages = {"index": "index.html", "publist": "publications.html"}

# If false, no module index is generated.
html_domain_indices = True

html_index = "index.html"

# If false, no index is generated.
html_use_index = True

html_use_modindex = True

# If true, the index is split into individual pages for each letter.
html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# does not have any affect. sphinx credit is located in footer.html
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
html_use_opensearch = "https://pygimli.org"

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = "gimlidoc"

# -- Options for LaTeX output

# The name of an image file (relative to this directory) to place at the top of
# the title page.
latex_logo = join(SPHINXDOC_PATH, "_static/gimli_logo.pdf")

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
# latex_use_parts = False

# If true, show page references after internal links.
latex_show_pagerefs = True

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True

extradir = path.abspath(join(SPHINXDOC_PATH, "_static"))  # .replace("\\", "/")

latex_elements = {
    # The paper size ("letterpaper" or "a4paper").
    "papersize": "a4paper",
    # The font size ("10pt", "11pt" or "12pt").
    "pointsize": "11pt",
    # Additional stuff for the LaTeX preamble.
    "preamble": "\\usepackage{amsfonts}\n\
    \\usepackage{amssymb}\n\
    \\usepackage{graphicx}\n\
    \\usepackage{amsmath}\n\
    \\usepackage{bm}\n\
    \\usepackage{pslatex}\n\
    \\graphicspath{{"
    + SPHINXDOC_PATH
    + "}}",
}

if sphinx.__version__.startswith("1.3"):
    latex_elements["preamble"] += "\\RequirePackage{fixltx2e}\n\
     \\MakeRobust\\DUspan\n"

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual])
latex_documents = [
    (
        "documentation",
        "gimli.tex",
        "GIMLi Documentation",
        r"Carsten Rücker, Thomas Günther \& Florian Wagner",
        "manual",
    )
]

try:
    pngmath_latex_preamble  # check whether this is already defined
except NameError:
    pngmath_latex_preamble = ""

pngmath_latex_preamble = "\
        \\usepackage{amsfonts}\n\
        \\usepackage{amssymb}\n\
        \\usepackage{amsmath}\n\
        \\usepackage{bm}\n\
        \\usepackage{pslatex}\n"

# latex.add_macro("\\newcommand{\\arr}[1]{\\ensuremath{\\textbf{#1}}}")

_mathpng_tempdir = DOC_BUILD_DIR + "mathtmp"

staticpath = os.path.abspath(join(SPHINXDOC_PATH, "_static"))
latex_additional_macros = open(join(staticpath, "mylatex-commands.sty"))

mathjax_latex_preamble = ""

LatexCommandTranslator = {}
LatexCommandTranslator["mathDictionary"] = {}
LatexCommandTranslator["commandDictionary"] = {}
LatexCommandTranslator["command1Dictionary"] = {}

trans = LatexCommandTranslator

for macro in latex_additional_macros:
    # used when building latex and pdf versions
    latex_elements["preamble"] += macro + "\n"
    # used when building html version
    pngmath_latex_preamble += macro + "\n"

    mathOperator = re.search("\\\\DeclareMathOperator{\\\\([A-Za-z]*)}{(.*)}", macro)
    if mathOperator:
        LatexCommandTranslator["mathDictionary"][mathOperator.group(1)] = (
            mathOperator.group(2).replace("\\", "\\\\")
        )

    newCommand = re.search("\\\\newcommand{\\\\([A-Za-z]*)}{(.*)}", macro)
    if newCommand:
        LatexCommandTranslator["commandDictionary"][newCommand.group(1)] = (
            newCommand.group(2).replace("\\", "\\\\")
        )
    newCommand = re.search("\\\\renewcommand{\\\\([A-Za-z]*)}{(.*)}", macro)
    if newCommand:
        LatexCommandTranslator["commandDictionary"][newCommand.group(1)] = (
            newCommand.group(2).replace("\\", "\\\\")
        )

    newCommand = re.search("\\\\newcommand{\\\\([A-Za-z]*)}\\[1\\]{(.*)}", macro)
    if newCommand:
        LatexCommandTranslator["command1Dictionary"][newCommand.group(1)] = (
            newCommand.group(2).replace("\\", "\\\\")
        )
    newCommand = re.search("\\\\renewcommand{\\\\([A-Za-z]*)}\\[1\\]{(.*)}", macro)
    if newCommand:
        LatexCommandTranslator["command1Dictionary"][newCommand.group(1)] = (
            newCommand.group(2).replace("\\", "\\\\")
        )

latex_additional_macros.close()

plot2rst_commandTranslator = LatexCommandTranslator

# -- Options for manual page output -------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [("index", "GIMLi", "GIMLi Documentation", ["GIMLi Group"], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False

# -- Options for Texinfo output -----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        "index",
        "GIMLi",
        "GIMLi Documentation",
        "GIMLi Group",
        "GIMLi",
        "Geophysical Inversion and Modelling Library",
        "Miscellaneous",
    ),
]

# Documents to append as an appendix to all manuals.
# texinfo_appendices = []

# If false, no module index is generated.
# texinfo_domain_indices = True

# How to display URL addresses: "footnote", "no", or "inline".
texinfo_show_urls = "footnote"

# -- Options for pybtex output ------------------------------------------------
# load our plugins for manual bibstyle

# temporary disable due to python3 pybtex quirks
for dist in pkg_resources.find_distributions(
    SPHINXDOC_PATH + "/_templates/pybtex_plugins/"
):
    pkg_resources.working_set.add(dist)

# End pybtex stuff

# -- Options for doxylink -----------------------------------------------------
doxylink = {
    "gimliapi": (join(DOXY_BUILD_DIR, "gimli.tag"), "https://www.pygimli.org/gimliapi")
}

# Create HTML table
from bib2html import write_html

publications = write_html()

################################################################################################
# Extra call to create small gallery of all already made tutorials and examples in the sidebar.
################################################################################################
make_gallery(os.path.abspath(SPHINXDOC_PATH), os.path.abspath(DOC_BUILD_DIR))

# Add carousel to start page
from paper_carousel import showcase

random.shuffle(showcase)  # mix it up
html_context = {
    "showcase": showcase,
    "publications": publications,
    "default_mode": "light",
}

srclink_project = "https://github.com/gimli-org/gimli"
srclink_src_path = "doc/"
srclink_branch = "dev"

################################################################################
# -- Options for Napoleon
################################################################################
# New docstring parsing using Napoleon instead of numpydoc
# The monkeypatch detects draw or show commands
from sphinx.ext.napoleon import NumpyDocstring

def monkeypatch(self, section: str, use_admonition: bool):
    lines = self._strip_empty(self._consume_to_next_section())
    lines = self._dedent(lines)
    all_lines = " ".join(lines)
    if (("show" in all_lines or "draw" in all_lines)
        and "Example" in section
        and ".. plot::" not in all_lines
       ):
        header = [f".. rubric:: {section}"]
        header.append(".. plot::\n\n")
        lines = self._indent(lines, 3)
    else:
        header = [f".. rubric:: {section}"]

    if lines:
        return header + [""] + lines + [""]
    else:
        return header + [""]

NumpyDocstring._parse_generic_section = monkeypatch

napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True

################################################################################
# -- Options for Bibtex
################################################################################

bibtex_bibfiles = ["gimliuses.bib", "libgimli.bib", "references.bib"]
bibtex_reference_style = "author_year"

################################################################################
# -- Options for breath
################################################################################

# breathe_projects = {"gimli":
#                         os.path.abspath(join(DOC_BUILD_DIR, "../../doxygen/xml")),
#                     }
# 
# breathe_default_project = "gimli"
# Settings for myst
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "substitution",
]
myst_dmath_allow_labels=True
nb_execution_excludepatterns = ["*Untitled*", "_examples_auto/**/*", "_tutorials_auto/**/*"]
