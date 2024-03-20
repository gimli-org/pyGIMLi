---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"editable": true, "slideshow": {"slide_type": ""}}

# User guide
```{contents} Table of Contents
:depth: 3
```

+++ {"editable": true, "slideshow": {"slide_type": ""}}

## Citations

+++ {"editable": true, "slideshow": {"slide_type": ""}}

{cite}`Seigel1959` did something great

+++ {"editable": true, "slideshow": {"slide_type": ""}}

## Evaluate code

```{code-cell} ipython3
---
editable: true
slideshow:
  slide_type: ''
---
import pygimli as pg
mesh = pg.createGrid(20,5)
data = pg.x(mesh)
pg.show(mesh, data)
```

##  Some Markdown features

### Math
Since Pythagoras, we know that {math}`a^2 + b^2 = c^2`.

$$
(a + b)^2  &=  (a + b)(a + b) \\
           &=  a^2 + 2ab + b^2
$$ (mymath2)

The equation {eq}`mymath2` is also a quadratic equation.

Some **text**!

:::{admonition} Here's my title
:class: tip

Here's my admonition content.
:::

### Tables
:::{table} Table caption
:widths: auto
:align: center

| foo | bar |
| --- | --- |
| baz | bim |
:::

### Typography
**strong**, _emphasis_, `literal text`, \*escaped symbols\*

### Footnotes
A longer footnote definition.[^mylongdef]

[^mylongdef]: This is the _**footnote definition**_.

    That continues for all indented lines

    - even other block elements

    Plus any preceding unindented lines,
that are not separated by a blank line

This is not part of the footnote.

### Cards
:::{card} Card Title
Header
^^^
Card content

+++

Footer
:::

### Tabs
::::{tab-set}

:::{tab-item} Label1
Content 1
:::

:::{tab-item} Label2
Content 2
:::

::::
