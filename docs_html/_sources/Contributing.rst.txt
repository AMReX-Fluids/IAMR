.. _Chap:Contributing:

Contributing
===============

Development Model
-----------------

Development generally follows the following ideas:

* New features are merged into to the ``development`` branch using
  Pull Requests (PRs).

  Nightly regression testing is used to ensure that no answers
  change (or if they do, that the changes were expected).

* Bug fixes, questions and contributions of new features are welcome!

 * Bugs should be reported through GitHub Issues.
 * We suggest asking questions through GitHub Discussions.
 * All contributions should be done via pull requests.
   A pull request should be generated from your fork of
   IAMR and target the ``development`` branch. See below for
   details on how this process works.

   In general we squash commits upon merge to have a clean history.
   *Please ensure that your PR title and first post are descriptive,
   since these will be used for a squashed commit message.*

   Please note the following:
     If you choose to make contributions to the code
     then you hereby grant a non-exclusive, royalty-free perpetual license
     to install, use, modify, prepare derivative works,
     incorporate into other computer software,
     distribute, and sublicense such enhancements or derivative works
     thereof, in binary and source code form.

* On the first workday of each month, we make a tagged release.  The merge window into
  ``development`` is closed a few days before the release day.  While the merge window is closed,
  only bug fixes should be merged into ``development``.  Once the release is done, the merge window
  reopens.

Git workflow
------------

IAMR uses `git <https://git-scm.com>`_ for version control. If you
are new to git, you can follow one of these tutorials:

- `Learn git with bitbucket <https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud>`_
- `git - the simple guide <http://rogerdudler.github.io/git-guide/>`_


The basic workflow is:

- Fork the main repo (or update it if you already created it).
- Implement your changes on a new branch ``<branch_name>`` on
  your fork. Don't forget to push.
- Create a Pull Request from branch ``<branch_name>`` on your fork to branch
  ``development`` on the main IAMR repository.


.. _sec:forking:

Make your own fork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
First, let us setup your local git repo. To make your own fork of the main
(`upstream`) repository, press the fork button on the `IAMR Github page <https://github.com/IAMR-Codes/IAMR>`_.

Then, clone your fork on your local computer. If you plan on doing a lot of IAMR development,
we recommend configuring your clone to use ssh access so you won't have to enter your Github
password every time, which you can do using these commands:

.. code:: shell

  git clone --branch development git@github.com:<myGithubUsername>/IAMR.git

  # Then, navigate into your repo, add a new remote for the main IAMR repo, and fetch it:  
  cd IAMR
  git remote add upstream https://github.com/IAMR-Codes/IAMR
  git remote set-url --push upstream git@github.com:<myGithubUsername>/IAMR.git
  git fetch upstream

  # We recommend setting your development branch to track the upstream one instead of your fork:
  git branch -u upstream/development


For instructions on setting up SSH access to your Github account on a new machine, see `here <https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.

If you instead prefer to use HTTPS authentication, configure your local clone as follows:

.. code:: shell

  git clone --branch development https://github.com/<myGithubUsername>/IAMR.git

  # Navigate into your repo, add a new remote for the main IAMR repo, and fetch it
  cd IAMR
  git remote add upstream https://github.com/IAMR-Codes/IAMR
  git remote set-url --push upstream https://github.com/<myGithubUsername>/IAMR.git
  git fetch upstream

  # We recommend setting your development branch to track the upstream one instead of your fork:
  git branch -u upstream/development
  
Note: you do not have to re-do the setup above every time.
Instead, in the future, to update the ``development`` branch on your fork
(i.e. sync your local ``development`` branch with that of the main repo) use

.. code:: shell

  git checkout development
  git pull
  
Now you are free to play with your fork (for additional information, you can visit the
`Github fork help page <https://help.github.com/en/articles/fork-a-repo>`_.



Making code changes
^^^^^^^^^^^^^^^^^^^^^^^^^^

We do `not` reccomend making any changes within the ``development`` branch. Doing so
will make your local ``development`` branch diverge from the matching branch in the
main IAMR repository.
Instead, we reccomend creating a new branch based on the ``development`` branch.
To proceed, make sure you are on the ``development`` branch within the IAMR directory.
Then create branch ``<branch_name>`` (the branch name should reflect the piece
of code you want to add, like ``high_order_interpolation``) with

.. code:: shell

  git checkout development
  git checkout -b <branch_name>

Do the coding you want.
Add the files you worked on to the git staging area with

.. code:: shell

  git add <file_I_created> <and_file_I_modified>

Periodically commit your changes with

.. code:: shell

  git commit -m "This is a 50-char description to explain my work"


The commit message (between quotation marks) is super important in order to
follow the developments and identify bugs.

For the moment, commits are on your local repo only. You can push them to
your fork with

.. code:: shell

  git push -u origin <branch_name>

If you want to synchronize your branch with the main IAMR ``development`` branch (this is useful
when ``development`` is being modified while you are working on
``<branch_name>``), you can use

.. code:: shell

  git merge upstream/development
  
and fix any conflicts that may occur.

Again, do not merge your branch into your local ``development`` branch,
because it will make your local ``development`` branch diverge from the
matching branch in the main repository. Instead, create a Pull Request (PR).

Submitting a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^^^

A Pull Request is the way to efficiently visualize the changes you made
and to propose your new feature/improvement/fix to the IAMR project.
Right after you push changes, a banner should appear on the Github page of
your fork, with your ``<branch_name>``.

- Click on the ``compare & pull request`` button to prepare your PR.
- It is time to communicate your changes: write a title and a description for
  your PR. People who review your PR are happy to know
  * what feature/fix you propose, and why
  * how you made it (created a new class than inherits from...)
  * and anything relevant to your PR (performance tests, images, *etc.*)
- Press ``Create pull request``. Now you can navigate through your PR, which
  highlights the changes you made.

Please DO NOT write large Pull Requests, as they are very difficult and
time-consuming to review. As much as possible, split them into small,
targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.
If you want to implement a large feature, write helper functionality first, test it and submit those as a first pull request.
If you want to implement a feature and are not too sure how to split it,
just open a discussion about your plans and ping other IAMR developers on it to chime in.

Even before your work is ready to merge, it can be convenient to create a PR
(so you can use Github tools to visualize your changes). In this case, please
make a ``draft`` PR using the drop-down menu next to the ``Create pull request`` button.

Once your pull request is made, we will review and potentially merge it.
We recommend always creating a new branch for each pull request, as per the above instructions.
Once your pull request is merged, you can delete your PR branch with

.. code:: shell

  git branch -D <branch_name> # Delete the local copy
  git push origin --delete <branch_name> # Delete the remote one on your fork

Generally speaking, you want to follow the following rules.

  * Do not merge your branch for PR into your local ``development`` branch that tracks IAMR
    ``development`` branch.  Otherwise your local ``development`` branch will diverge from IAMR
    ``development`` branch.

  * Do not commit in your ``development`` branch that tracks IAMR ``development`` branch.

  * Always create a new branch based off ``development`` branch for each pull request, unless you are
    going to use git to fix it later.

If you have accidentally committed in ``development`` branch, you can fix it as follows,

.. code:: shell

  git checkout -b new_branch
  git checkout development
  git reset HEAD~2  # Here 2 is the number of commits you have accidentally committed in development
  git checkout .

After this, the local ``development`` should be in sync with IAMR ``development`` and your recent
commits have been saved in ``new_branch`` branch.

If for some reason your PR branch has diverged from IAMR, you can try to fix it as follows.  Before
you try it, you should back up your code in case things go wrong.

.. code:: shell

  git fetch upstream   # assuming upstream is the remote name for the official IAMR repo
  git checkout -b xxx upstream/development  # replace xxx with whatever name you like
  git branch -D development
  git checkout -b development upstream/development
  git checkout xxx
  git merge yyy  # here yyy is your PR branch with unclean history
  git rebase -i upstream/development

You will see something like below in your editor,

.. code:: shell

  pick 7451d9d commit message a
  pick c4c2459 commit message b
  pick 6fj3g90 commit message c

This now requires a bit of knowledge on what those commits are, which commits have been merged,
which commits are actually new.  However, you should only see your commits.  So it should be
easy to figure out which commits have already been merged.  Assuming the first two commits have been
merged, you can drop them by replace ``pick`` with ``drop``,

.. code:: shell

  drop 7451d9d commit message a
  drop c4c2459 commit message b
  pick 6fj3g90 commit message c

After saving and then exiting the editor, ``git log`` should show a clean history based on top of
``development`` branch.  You can also do ``git diff yyy..xxx`` to make sure nothing new was dropped.  If
all goes well, you can submit a PR using ``xxx`` branch.
Don't worry, if something goes wrong during the rebase, you an always ``git rebase --abort`` and start over.


IAMR Coding Style Guide
-----------------------

Code Guidelines
^^^^^^^^^^^^^^^^^

IAMR developers should adhere to the following coding guidelines:

* Indentations should use 4 spaces, not tabs.
* Use curly braces for single statement blocks. For example:

  .. code::
       
     for (int n=0; n<10; ++n) {
         Print() << "Like this!";
     }
     
     for (int n=0; n<10; ++n) { Print() << "Like this!"; }

  but not
  
  .. code::

     for (int n=0; n<10; ++n) Print() << "Not like this.";
     
     for (int n=0; n<10; ++n)
        Print() << "Not like this.";

* Add a space after the function name and before the
  parenthesis of the parameter list (but
  not when simply calling the function). For example:

  .. code::
     
     void CorrectFunctionDec (int input)

  Not

  .. code::

     void IncorrectFunctionDec(int input)

  This makes it easy to find where functions are defined with grep.
* Member variables should be prefixed with ``m_``. For example:

  .. code::
     
     IAMR::Real m_variable;

These guidelines should be adhered to in new contributions to IAMR, but
please refrain from making stylistic changes to unrelated sections of code in your PRs.

API Documentation Using Doxygen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Doxygen documentation is designed for advanced user-developers. It aims
to maximize the efficiency of a search-and-find style of locating information.
Doxygen style comment blocks should proceed the namespace, class, function, etc.
to be documented where appropriate. For example:

.. code:: 

  /**
   * \brief A one line description.
   *
   * \param[in] variable A short description of the variable.
   * \param[inout] data The value of data is read and changed.
   *
   * A longer description can be included here.
   */

    void MyFunction (int variable, MultiFab& data){
    ...

Additional information regarding Doxygen comment formatting can be found
in the `Doxygen Manual <https://www.doxygen.nl/manual/>`_.
