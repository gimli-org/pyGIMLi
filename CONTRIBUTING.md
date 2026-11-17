Contributing
============

As an open-source project, **pyGIMLi always welcomes contributions from
the community**. Here, we offer guidance for 3 different ways of
contributing with increasing levels of required coding proficiency.

A. Submit a bug report
----------------------

If you expericence issues with pyGIMLi or miss a certain feature, please
[open a new issue on GitHub](https://github.com/gimli-org/gimli/issues).
To do so, you need to [create a GitHub
account](https://github.com/join).

B. Send us your example
-----------------------

We are constantly looking for interesting usage examples of pyGIMLi. If
you have used the package and would like to contribute your work to the
chapt:examples, please send your script to <mail@pygimli.org>. Make sure
that the individual steps in your Python script are documented according
to the [sphinx-gallery
syntax](http://sphinx-gallery.readthedocs.io/en/latest/tutorials/plot_notebook.html).

C. Contribute to the code
-------------------------

> **note**
>
> To avoid redundant work, please [contact us](mailto:mail@pygimli.org)
> before you start working on a non-trivial feature.

The preferred way to contribute to the pygimli code base is via a pull
request (PR) on GitHub. The general concept is explained
[here](https://guides.github.com/introduction/flow) and involves the
following steps:

### 1. Fork the repository

If you are a first-time contributor, you need [a GitHub
account](https://github.com/join) and your own copy ("fork") of the
code. To do so, go to <https://github.com/gimli-org/gimli> and click the
"Fork button" in the upper right corner. This will create an identical
copy of the complete code base under your username on the GitHub server.
Clone this repository to your local disk:

``` {.sourceCode .bash}
git clone https://github.com/YOUR_USERNAME/gimli
```

After that you can install the software as usual (see sec:install).

### 2. Create a feature branch

Go to the source folder and create a feature branch to hold your
changes. It is advisable to give it a sensible name such as
`adaptive_meshes`.

``` {.sourceCode .bash}
cd gimli
git checkout -b adaptive_meshes
```

### 3. Start making your changes

Go nuts! Add and modify files and regularly commit your changes with
meaningful commit messages. Remember that you are working in your own
personal copy and in case you break something, you can always go back.
While coding, we encourage you to follow a few sec:coding\_guidelines.

``` {.sourceCode .bash}
git add new_file1 new_file2 modified_file1
git commit -m "Implemented adaptive meshes."
```

### 4. Test your code

Make sure that everything works as expected. New functions should always
contain a docstring with a test:

``` {.sourceCode .python}
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

When you run `pg.test()` the docstring test will be evaluated. See also
the section on sec:testing.

### 5. Submit a pull request

Once you implemented a functioning new feature, make sure your GitHub
repository contains all your commits:

``` {.sourceCode .bash}
git push origin adaptive_meshes
```

After pushing, you can go to GitHub and you will see a green PR button.
Describe your changes in more detail. Once reviewed by the core
developers, your PR will be merged to the main repository.

### 6. Updating your work with changes from upstream

While you work on your forked repository, sometimes changes are committed to the
main repository (usually called **upstream**).
You do NOT need to delete your forked repository and refork to apply these
changes to your own fork.
Follow the procedure described `here
<https://help.github.com/en/articles/syncing-a-fork>`__

Required only once, add the main repository as a remote to your local (cloned)
git repository:

   git remote add upstream https://github.com/gimli-org/gimli.git

Development work is usually based on the *dev* branch, and before creating a
new feature branch, make sure to update this branch on your forked repository:

.. code:: bash

	git checkout -b dev
	git merge upstream/dev
	# the following will push updates to your personal fork
	git push

Only then create your feature branch

.. code:: bash

   git checkout -b new_feature
