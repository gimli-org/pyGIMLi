---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Cheatsheet

+++

```{contents} Table of Content
:depth: 3
```

+++

## Citations

+++

{cite}`Seigel1959` did something great

+++

## Evaluate code

```{code-cell} ipython3
import pygimli as pg
mesh = pg.createGrid(20,5)
data = pg.x(mesh)
pg.show(mesh, data)
```

## Some Markdown features

### Math

Since Pythagoras, we know that $a^2 + b^2 = c^2$

$$
(a + b)^2  =  (a + b)(a + b) \\
           =  a^2 + 2ab + b^2
$$ (mymath2)

The equation {eq}"mymath2" is also a quadratic equation.


## Custom latex macros

Some custom latex macros inspired by [physics latex package](https://ctan.org/):

| Macro       | Output         |
|-------------|----------------|
| \order{h^{2}} | $\order{h^{2}}$ |
| \vb{v}        | $\vb{v}$        |
| \grad{u}    | $\grad{u}$    |
| \div{\vb{v}}  | $\div{\vb{v}}$  |
| \curl{\vb{v}}        | $\curl{\vb{v}}$ |
| \laplacian{u}   | $\laplacian{u}$   |
| \sin(x)        | $\sin(x)$        |
| \dd x         | $\dd x$         |


Some **text-like stuff**!

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
