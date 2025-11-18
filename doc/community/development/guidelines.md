(sec-coding-guidelines)=

# Developer guidelines

## Contributing to the code

:::{note}
These steps are only relevant if you have write privileges to the main
GitHub repository. The general contribution guidelines can be found in
{ref}`sec-contributing`.
:::

1. Clone the GIMLi repository:

   > ```bash
   > git clone https://github.com/gimli-org/gimli gimli-src && cd gimli-src
   > ```

You should be on the master branch by default. Check by calling *git status*.
GIMLi encourages developers to keep the *master branch* clean. So please use the
default development branch called *dev*. Save [this
script](https://gist.github.com/florian-wagner/0be0b09bbabec6791ed3) in
.git/hooks for automatic warnings in case you are about to push to the master
branch.

2. Checkout the dev branch:

   > ```bash
   > git checkout dev
   > ```

3. Make changes in your favorite editor.

4. Add new files to the staging area:

   > ```bash
   > git add new_file1 new_file2
   > ```

5. Make a commit with a *meaningful* message:

   > ```bash
   > git commit -a -m "ADD: important new feature."
   > ```
   >
   > Make use of the prefixes in the table below to categorize commits.

6. Pull the latest developments from GitHub using automatic rebase:

> For more info see: <http://kernowsoul.com/blog/2012/06/20/4-ways-to-avoid-merge-commits-in-git/>
>
> > ```bash
> > git pull --rebase
> > ```

You can set the rebase behavior on default with:

> ```bash
> git config --global branch.autosetuprebase always
> ```

7. Push to the origin development branch:

   > ```bash
   > git push origin dev
   > ```

Note that if you are on the `dev` branch, a `git push` should suffice, but the
command above is more explicit about the branch which should be pushed.

8. Show available tags and branches:

   > ```bash
   > git tag
   > git branch
   > ```

9. If you plan bigger changes or developments create an own branch for it:

   > ```bash
   > git branch my_new_test_branch
   > git checkout my_new_test_branch
   > ```

Work on your code in branch `my_new_test_branch` and test it.
If you are happy with your results merge them back
into the dev branch and delete the branch.

> ```bash
> git checkout dev
> git merge my_new_test_branch
> git branch -d my_new_test_branch
> ```

## Branches

| Branch name | Description             | Purpose                                                                                                                                                                                                                               | CI (Jenkins) Current behavior                                                                                                                      | CI (Jenkins) Wishlist                                                                                                        | Rules                                                                                                                 |
| ----------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| dev         | Main development branch | This is where new versions are developed and pull requests are merged to.                                                                                                                                                             | Build and test, on success generate <http://dev.pygimli.org> and merge to master                                                                   | Test pull requests, build nightly conda packages for test channel                                                            | Run pg.test() and make doc before pushing. Tag related issues in commit message.                                      |
| release     | Latest release          | This is where releases are "staged". This usually means creating a git tag and manually merging dev into release. Hot fixes and website typos can be directly committed here.                                                         | Start in empty workspace, build, test and create documentation after each push. On success, merge into master and update <http://www.pygimli.org>. | Test "backward compatibility" (e.g., run example scripts from last release with this version). Test on Mac and Windows, too. | Make sure the tag is annotated and the version string is following the principles described below.                    |
| master      | Latest tested dev       | Make sure that if people checkout the repository, they always have a working version.                                                                                                                                                 | (represented by <http://dev.pygimli.org>)                                                                                                          | Build pgcore (if necessary) and pygimli conda packages for release.                                                          | Never push anything to master!                                                                                        |
| -           | Feature branches        | Larger endeavors and major restructuring should happen in dedicated feature branches (or forks), which are eventually merged to dev. This can also be useful if you want to give write access to others to jointly work on a feature. | Automatic testing can be requested (<mailto:florian@pygimli.org>).                                                                                 |                                                                                                                              | Start feature branch from dev. Inform other developers about your develpment (to avoid conflicts and redundant work). |

## Commit messages

Please prepend your commit messages with one of the following prefixes depending on the nature of the commit:

| **Prefix** | **Meaning**                                                          |
| ---------- | -------------------------------------------------------------------- |
| **ADD:**   | Addition of new functionality                                        |
| **API:**   | API change (argument orders or renames)                              |
| **BLD:**   | Changes to pyGIMLi's build pipeline                                  |
| **CI:**    | Continuous integration (i.e., Jenkins related)                       |
| **CLN:**   | Clean up, refactoring or typo correction                             |
| **DEF:**   | Change of default behavior                                           |
| **DEP:**   | Deprecation                                                          |
| **DOC:**   | Improve documentation (i.e., docstrings or website)                  |
| **ENH:**   | Enhancement, e.g. due to more flexibility by new (keyword) arguments |
| **FIX:**   | Fixing a bug                                                         |
| **REL:**   | Release (setting tags, updating version strings, etc.)               |
| **STY:**   | Coding style improvement (PEP8, PEP257)                              |
| **TST:**   | New or improved test                                                 |

:::{note}
The tags **API, ENH, ADD, FIX** are relevant for creating the changelog later.
:::

## Version numbering

Following the PEP conventions, pygimli holds a `__version__` string.

```python
import pygimli as pg
print(pg.__version__)
v1.4.1-0-g76a4848-modified
```

This string consists of:

> \<*Last git tag*>-\<*Number of commits since last tag*>-\<*last commit hash*>

The string has a *-modified* suffix if pygimli has uncommitted changes.

To produce a new version, type:

```bash
git tag -a -m "First official release" "v1.0.x" # tags last commit as v1.0.x
git push --tags # pushes tags to GitHub
```

(sec-testing)=

To see the commits since the last tag/release, you can use:

```bash
git log v1.5.0...HEAD --oneline
```

Or to see commits between specific versions:

```bash
git log v1.5.3...v1.5.4 --oneline
```

Alternatively, this information can also be obtained [via GitHub](http://github.com/gimli-org/gimli/compare/v1.4.0...v1.4.1).

## Testing

Run specific API examples from shell:

```bash
python -c "import pygimli as pg; pg.test(pg.meshtools.createCircle, show=True)"
```

Run a specific test from shell:

```bash
python -c "import pygimli; from pygimli.physics.petro.resistivity import *; test_Archie()"
```

Run all tests:

```bash
python -c "import pygimli; pygimli.test(show=True)"
```

Run pylint from shell to check code:

```bash
pylint --rcfile $GIMLIROOT/.pylintrc file.py
```

Run prospector to check code like landscape.io do:

```bash
prospector --profile=$GIMLIROOT/.prospector.yml file.py
```

Read API documentation from shell:

```bash
python -c "import pygimli as pg; help(pg.test)"
```

### Adding an example to the paper carousel

1. Put a catchy figure in *doc/\_static/usecases*. It should be 750X380 px. You can use this command:

```bash
mkdir -p converted
mogrify -path converted -format jpg -resize "750x380^" -gravity North -crop 750x380+0+0 +repage -quality 95 *.jpg
```

2. Add a new dictionary with additional information into *doc/paper_carousel.py*. The format should be like:

```python
dict(img="jordi2018.jpg",
     title="Improved geophysical images through geostatistical regularization",
     subtitle="Jordi et al. (2018), Geophysical Journal International",
     link="https://doi.org/10.1093/gji/ggy055")
```

(sec-coding-rules)=

## Coding Guidelines

- Generally, we use PEP 8 <https://www.python.org/dev/peps/pep-0008/>?
- All names should be literally and in CamelShape style.
- Classes starts with upper-case letters, members and methods start with lower-case letters.
- All class members (self.member) need to be initialized in the init.
- (ugly²) Do not use multi initialize in one line, e.g., `a, b, c = 0, 0, 0`
- check for data types with `if isinstance(var, type):` instead `if type(var) == type:`

(sec-coding-guidelines-1)=

Use pylint or prospector to improve code quality.

We use: (exceptions in .landscape.yml and .pylintrc)

- pep8
- pep257
- pylint
- pyflakes

Add two blank lines between classes and methods but never use two blank lines otherwise.

### Behavior by name for global functions

```python
createFOO(...)
    """Return an instance of FOO."""
```

```python
showFOO(Bar, ...)
    """Open a window or optionally use a given xxes to show Bar as Foo."""
    return ax, cbar
```

```python
drawFOO(ax, Bar...)
    """Always needs an axes object ax first argument to draw Bar as Foo."""
    return graphics_object
```

```python
readFOO(fileName, *args):
    """Read object from disc."""
    return obj
```

```python
importFOO(fileName, obj, *args):
    """Import object from disc into an existing object."""
    return obj
```

```python
exportFOO(obj, fileName):
    """Export object to disc in foreign (FOOs) format."""
    return true
```

```python
convertFOO(fooObj):
    """Convert Foo obj into pyGIMLi Obj."""
    return gimliObj
```

### API Documentation and doctests

Use the following documentation syntax or see:
<https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html>

```python
def foo(arg1:type, arg2:type=Default, *args, **kwargs) -> ReturnType:
"""Short description, i.e., one line to explain what foo does. [PERIOD_AT_END]
  [ONE BLANKLINE]
  Explain a little more verbose was foo does. Use references :cite:`Archie1942`
  Use links to pygimli api :gimliapi:`GIMLI::Mesh`, `:py:mod:pygimli.manager`
  for modules, `:py:func:pygimli.solver.solveFiniteElements` for functions

  Type hints can be used for arguments and return values.

  Use math.

  .. math :
      a + \sigma \cdot \rho

  Explain all parameters.

  Arguments
  ---------
  arg1: type=Default | use links to :gimliapi:`GIMLI::Mesh`
      Describe arg1.
  arg2: type
      Describe arg2.

  Keyword Args
  ------------
  args: type
      Description.

  Attributes
  ----------
  For members

  Returns
  -------
  type:
      Describe the return value. If it is a class instance, use links to
      :gimliapi:`GIMLI::Mesh` or `:py:class:pygimli.core.Mesh`.

  Examples
  --------
  >>> import foo
  >>> ...
  >>> ...
  expected output  # for doctest to pass

  See Also
  --------
      average: Weighted average,
      e.g., Link to tutorials :ref:`tut:Modelling_BC` assuming
      an appropriate label in the tutorial.

  References
  ----------
  If not already in the global bib file.
  """
```
