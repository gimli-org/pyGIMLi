(sec-contributing)=

# How to contribute

As an open-source project, pyGIMLi always welcomes contributions from the
community. Here, we offer guidance for four different ways of contributing with
increasing levels of required coding proficiency (A-D).

## A. Submit a bug report

If you experience issues with pyGIMLi or miss a certain feature, please [open a
new issue on GitHub](https://github.com/gimli-org/gimli/issues). To do so,
you need to [create a GitHub account](https://github.com/join).

## B. Send us your example

We are constantly looking for interesting usage examples of pyGIMLi. If you have
used the package and would like to contribute your work to the {ref}`chapt:examples`, you can send us your example.
Make sure that your Python script is documented according to the [sphinx-gallery syntax](https://sphinx-gallery.github.io/stable/tutorials/plot_parse.html#sphx-glr-tutorials-plot-parse-py).

## C. Quickly edit the documentation

If you would like to edit a typo or improve the documentation on our website, you can submit fixes to the
documentation page by clicking on "Improve this page" under the quick links menu on every page in the right sidebar.
This will direct you to the Github page to make the edits directly on the ".rst" using Github (you will need a Github account, see step A).
Once you made your edits, fill out the fields under "Commit changes" (provide an overview of the changes you have made).
Note that you will make changes to a new branch under your fork since you will not have access to make changes in the existing branch.
Once you click "Commit changes", this automatically creates a pull request which we will review and then merge, if everything is good to go!

If you would like to contribute to the API or make changes offline to the files under `doc`, see following *Contribute to Code* and *Documentation* for more information.

## D. Contribute to code

:::{note}
Before you start working on a non-trivial feature, consider [contacting us](mailto:mail@pygimli.org)
or check the [current issues](https://github.com/gimli-org/gimli/issues) to avoid redundant work
:::

The preferred way to contribute to the pygimli code base is via a pull request
(PR) on GitHub. The general concept is explained [here](https://guides.github.com/introduction/flow) and involves the following steps:

### 1. Fork and clone the repository

If you are a first-time contributor, you need [a GitHub account](https://github.com/join) and your own copy ("fork") of the code. To do so, go
to <https://github.com/gimli-org/pyGIMLi> and click the "Fork button" in the upper right corner.
This will create a copy of the complete code base under your username on the GitHub server.
You can navigate to this repository and clone it to your local disk:

```bash
git clone https://github.com/YOUR_USERNAME/pyGIMLi
```

### 2. Prepare your environment

To install the neccesary Python requirements on your computer, we recommend the conda package manager in a lightweight installation such as Miniforge.
The cloned and forked repository contains an *environment.yml* file with the specification for all package requirements to build and test pyGIMLi.
This file should exist in the directory where you have cloned the repository. In order to create a separate environment where you will work,
run the following command in the same directory of the repository in your terminal:

```bash
conda env create # Creates the environment (only once)
conda activate pgdev # Activates it (everytime you want to work in it)
conda develop . # Installs your git-version of pygimli (only once)
python -c "import pygimli; pygimli.test()" # Make sure that everything works (also after adding new code)
```

You will need to do the last step (`conda activate pgdev`) everytime you start a terminal (or put it in your `.bashrc` startup file).
The prompt `(pgdev)` should appear before your terminal location.

### 3. Create a feature branch

Go to the source folder and create a feature branch to hold your changes. It is
advisable to give it a meaningful name, e.g. `adaptive_meshes`.

```bash
cd gimli
git checkout -b adaptive_meshes
```

### 4. Start making your changes

Go nuts! Add and modify files and regularly commit your changes with meaningful
commit messages. Remember that you are working in your own personal copy and in
case you break something, you can always go back. While coding, we encourage you
to follow a few {ref}`sec-coding-rules`.

```bash
git add new_file1 new_file2 modified_file1
git commit -m "ADD: Implemented adaptive meshes."
```

### 5. Test your code

Make sure that everything works as expected. New functions should always contain
a docstring with a test:

```python
def sum(a, b):
    """Return the sum of `a` and `b`.

    Examples
    --------
    >>> a = 1
    >>> b = 2
    >>> sum(a,b)
    3
    """
    return a + b
```

When you run `pg.test()` the docstring test will be evaluated. See also the
section on {ref}`sec-testing`.

### 6. Documentation

We use sphinx to build the web pages from these sources. To edit HTML file, navigate to the doc folder in the repository.
Once you're done making edits, run the following commands on your terminal.

```bash
make html
```

The process should build all of the documentation so it may take some time to complete. When it finishes, the last line should state:

```bash
Build finished. The HTML pages are in _build/html.
```

### 7. Submit a pull request

Once you implemented a functioning new feature, make sure your GitHub repository
contains all your commits:

```bash
git push origin adaptive_meshes
```

After pushing, you can go to GitHub and you will see a green PR button. Describe
your changes in more detail. Once reviewed by the core developers, your PR will
be merged to the main repository.

:::{note}
Please see {ref}`sec-coding-guidelines` for additional notes on coding rules, testing, etc.
:::
