.. _contribute:

==========
Contribute
==========

Requirements
============

1. make sure your local ``main`` branch is up-to date (``git switch main && git pull``)
2. install the package in developer mode (see :ref:`installation-dev`)
3. install the git hooks (``pre-commit install``)

Development guidelines
======================

* The **Git workflow** used by this project is the `Feature branch workflow <https://docs.gitlab.com/ee/gitlab-basics/feature_branch_workflow.html>`__
* The Python **Code style** follows `PEP8 <https://peps.python.org/pep-0008/>`__
* The **docstring style** follows `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`__
* **Static code analysis** is performed with `ruff <https://beta.ruff.rs/docs/>`__
* **Unit-testing** is performed using `pytest <https://docs.pytest.org/en/latest/>`__
* **don't re-invent the wheel**, always check if a dependency has what you need before implementing it yourself
* **make clear atomic commits**,

  - make `atomic commits <https://www.aleksandrhovhannisyan.com/blog/atomic-git-commits/>`
    with a clear descriptive message
    (`best practices <https://initialcommit.com/blog/git-commit-messages-best-practices>`_)
  - always test your local branch *before* pushing

pre-commit
==========

The project comes with a pre-commit configuration which will run some
checks each time you attempt to commit.

Some checks will be applied automatically by the associated tool,
others with be applied but not staged yet so you can review them
before trying to commit again.

Make sure you have installed the git hooks before starting to commit
(``pre-commit install``).

Documentation
=============

| This documentation is based on
  `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`__.
| It makes use of
  `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__
  as the markup language within the
  `Furo <https://github.com/pradyunsg/furo>`__ Sphinx theme.

The automatic deployment of the documentation happens only on the protected branch (``main``).

While working on your contribution you can **build your development version** locally,

* ``make -C docs clean`` if you want to reset from previous builds
* ``make -C docs html``
* open ``docs/build/html/index.html`` with your favorite browser

Please, make sure that the documentation builds without warnings or errors,
as the documentation build is set to fail otherwise (in particular, warnings
are treated as errors, also by the project's Sphinx configuration).

Docstrings
----------

Each API member **has** to be documented with a docstring.

A simple template for a function is the following,

.. code-block:: rst

    A one-line summary that does not use variable names or the function name.

    A few sentences giving an extended description.

    Parameters
    ----------
    x : `~show.only.the_class`
        Description of parameter `x`.
    y : `int`
        Description of parameter `y`.

    Returns
    -------
    result : `show.the.full.class.path`
        Description of ``result``.

Please, see the `numpydoc style guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_ for more details.

Changelog
---------

This project is supported by a changelog to let users know what has been fixed
or improved.

Whenever you add a new contribution, you should always open first a new issue to track it,
then a Pull Request (PR) that closes that issue.

Every PR needs to add at least one *news fragment*, a small *.rst* file under ``docs/source/changes``
containing a description of what was changed in a way that is usefull for the user.

The name of the file has to follow the template ``PULL_REQUEST_ID.CONTRIBUTION_TYPE.rst``,
where ``CONTRIBUTION_TYPE`` has to be one of the following:

- removed
- deprecated
- added
- changed
- fixed

For example: a bug fix was found and an issue was opened.
A PR with ID #5 was opened to fix the bug and it introduced the file ``5.fixed.rst``
which reads,

.. code-block::

    Fixed a terrible bug which caused this function to output this nonsense result.

.. _contribute_auth:

Authorship
==========

The online git host service will take records of your contribution to the repository,

In case you commit with more than one email, please update the ``.mailmap``
file at the root of the project
(see `this page <https://www.git-scm.com/docs/gitmailmap>`__ on how to do it).

If you contribute to this code, please update also the ``CODEOWNERS`` file.
