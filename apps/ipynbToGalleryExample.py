#!/usr/bin/env python
r"""Convert jupyter notebook to sphinx gallery notebook styled examples.

Modified from here:https://gist.github.com/chsasank/7218ca16f8d022e02a9c0deb94a310fe#file-ipynb_to_gallery-py

Usage: python ipynbToGalleryExample.py <notebook.ipynb>

Dependencies:
pypandoc: install using `pip install pypandoc`
"""
import pypandoc as pdoc
import json
from pathlib import Path


def fixRST(rst):
    """Convert ipynb flavour reST into flat reST for sphinx-gallery."""
    #rst = rst.replace('\n', '') # only for windows
    rst = rst.replace('\r', '\n')
    rst = rst.replace(r":raw-latex:`\begin{align}", "\n.. math::\n\t")
    rst = rst.replace(r"\end{align}`", "")

    # some spocky encoding problems with '
    rst = rst.replace('’', "'")

    # We prefer $h$ instead of :math:`h` which seems valid for ipynb
    while ":math:`" in rst:
        tkLen = len(":math:`")
        start = rst.find(":math:`")
        end = rst.find("`", start + tkLen)
        rst = rst[:start] + '$' + rst[start + tkLen:end] + '$' + rst[end + 1:]

    return rst


def convert_ipynb_to_gallery(fileName):
    """Write the cells of an ipynb file as sphinx-gallery python script."""
    outFileName = fileName.replace('.ipynb', '.py')
    print(f"Converting {fileName} -> {outFileName}")

    header = "#!/usr/bin/env python\n" +\
             "# -*- coding: utf-8 -*-\n"

    with Path(fileName).open(encoding="utf-8") as fi:
        nb_dict = json.load(fi)

    cells = nb_dict['cells']

    first = True
    for cell in cells:
        try:
            if cell['source'][0].startswith('%'):
                continue
        except (IndexError, KeyError):
            pass

        if first:
            if cell['cell_type'] != 'markdown':
                continue
            first = False

            md_source = ''.join(cell['source'])

            rst_source = pdoc.convert_text(md_source, 'rst', 'md')
            python_file = header + 'r"""\n' + fixRST(rst_source) + '\n"""'
        else:
            if cell['cell_type'] == 'markdown':
                md_source = ''.join(cell['source'])
                rst_source = pdoc.convert_text(md_source, 'rst', 'md')
                rst_source = fixRST(rst_source)

                commented_source = '\n'.join(['# ' + x for x in
                                              rst_source.split('\n')])

                python_file = python_file + '\n\n' + '%%%'+ '\n' + \
                    commented_source

            elif cell['cell_type'] == 'code':
                source = ''.join(cell['source'])
                python_file = python_file + '\n\n# %%%\n' + source

    python_file = python_file.replace("\n%", "\n# %")
    Path(outFileName).write_text(python_file, encoding="utf-8")


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: ' + sys.argv[0] + ' ipythonnotebook.ipynb')
    else:
        convert_ipynb_to_gallery(sys.argv[-1])
