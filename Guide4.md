# Undoing things


## Introduction

In a perfect world, you would contribute code to a project and that code would be bug free. Similarly, you would make changes to your local files and the changes will bring nothing but joy.  This is often not the case. You may introduce new "features" that inadvertently break other code, or you may ignore a comment in your code, e.g. "Seriously, don't change this", change the code and break things. Let's imagine you've already comitted your buggy changes and realize your mistake.  How can you use ``git`` to revert back to a point where you know code worked?

Git has many tools to your disposal, and the general ``undoing`` tools are ``revert`` and ``reset``.  Superficially, both give the same result, but work in different ways.  In this guide, we'll go over the differences and usage cases of the aforementioned commands.


## The scenario

You are a new developer to a project, AlphaStar. You contribute code and everything is going smoothly.  Changes are comitted by you, pushed to the repository and everyone rejoices.  You continue the development and add a new feature to the file ``analyzer.py``.  After a couple of commits, you realize that this new code breaks a lot of other parts of the project. You decide that the best course of action is to revert the changes back to the original. Because you have been diligent, you kept a good commit history explaining the changes you made to ``analyzer.py`` with every commit. Great.  So you have two choices, before you push your commits back to the repository.  You can:

	* Change everything back to what is was, without leaving a trace
	* Change everything back to what is was while keeping track of the changes you made

The first scenario is typically for a case when you were being silly and it is irrelevant to keep a copy of these changes.  The second scenario is for when you would like to show that ``including feature breaks the code`` and to keep track of parts of the code that could (for instance) be fragile to changes. Typically, when you decide to rewrite the past, you would stick to the first scenario.  For this, you would use ``git reset``.  

### Git ``reset``

Git reset has different usage scenarios.  The most trivial of which, is the unstaging of files. To refresh, before you commit your changes to a repository, you first have to add a file to the staging area.

File1 -> edit1 -> add -> commit1
File1 -> edit2 -> add File1 -> commit2


You can ``reset`` the file to the last commit version in the _staging_ area.  This means that the current file you were working on will not be altered in any way.  So if we did a reset after "add File1" before ``commit2``, the repository will consider ``File1`` to be the ``commit1`` version. In the shell:

``shell
#Edit File1
$ > git add File1
$ > git commit -m "commit1"
#Edit File1
$ > git add File1
$ > git status
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

        modified:   File1

$ > git reset 
$ > git status
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   initial
```

You can see that after resetting, git is still aware of the file with changes and notifies you that the file has not been staged for commit.  

What if you really just want to overwrite your current version of ``File1`` with the one from  the most recent commit, ``commit1``? For this, you can use the (very dangerous) ``git reset --hard``.  This will overwrite any file you've edited with the version from ``commit1``.  Be _*very*_ careful when you use this command.  It _cannot_ be undone.


```shell
$ > git reset --hard
$ > git status
nothing to commit, working directory clean
```

### Git ``revert``

Git ``revert`` is a curious command.  As stated before, it is used to "revert" changes you've made and add it as a new commit. In our simplest case, the most recent commit (``commit2``) can be reverted:

```shell
$ > git log --oneline
6c34756 commit2
8b1b5f3 commit1
61da43a Initial commit

$ > git revert HEAD
```

The ``HEAD`` points to the most recent commit. What you should note, is that ``revert`` does not revert _*to*_ the ``commit2`` point (6c34756), it reverts _*that*_ commit, so we're back at ``commit1`` (8b1b4f3). 

The tricky bit of revert is that if we reverted ``commit1`` from the start - 

```shell
$ > git revert 8b1b5f3
```

- the general assumption is that this reverts _all_ commits from ``commit1`` onwards.  This is _not_ the case.  It attempts to revert the changes brought about by ``commit1``.  This is a little tricky, so lets create a new repository. 

```shell

$ > mkdir testrevert
$ > cd testrevert
$ > touch myfile.txt
$ > git add myfile.txt
$ > git commit -m "Initial commit"
$ > git log --oneline
c796082 Initial commit
```

Open ``myfile.txt`` and add the following text to it and save the file.

```
Line 1
Line 2
Line 3
Line 4
```

Add the ``myfile.txt`` and commit the changes with the message: "Our initial myfile.txt". In this simplistic example, the numbers for some of the lines will be replaced by the corresponding words (1 -> one, 2 -> two etc).  Edit the file and change ``Line 1`` to ``Line one``.  Save, add ``myfile.txt`` to the repository and commit with the message "Edited Line 1".  Do the same for ``Line 3`` with the commit message similar to the last. Your ``git log`` should look something like this:

```shell
f8e297f Edited Line 3
0c03424 Edited Line 1
a2ee513 Our initial myfile.txt
c796082 Initial commit
```

Let's revert the changes of the commit ``Edited Line 1``.  Look a the corresponding hash for your repo, and run:

```shell
$ > cat myfile.txt
Line one
Line 2
Line three
Line 4

$ > git revert [your "Edit Line 1" commit hash]
#It will ask you to leave a message giving a reason for this merge.  Save and exit the text editor
$ > git log
a3ed063 Revert "Edited Line 1"
f8e297f Edited Line 3
0c03424 Edited Line 1
a2ee513 Our initial myfile.txt
c796082 Initial commit

$ > cat myfile.txt
Line 1
Line 2
Line three
Line 4
```

You can do the same for the commit where you changed ``Line 1``.  You can even be silly and revert the revert. Git only sees commits as commits, so:

```shell
$ > git revert a3ed063
$ > git log --oneline
33c4185 Revert "Revert "Edited Line 1""
a3ed063 Revert "Edited Line 1"
f8e297f Edited Line 3
0c03424 Edited Line 1
a2ee513 Our initial myfile.txt
c796082 Initial commit
```

You will notice that the commit for ``Line 3`` has not been reverted.  This is because ``revert`` reverts changes of a specific commit and does not revert changes by the other commits.  Let's roll back a bit and revert all the changes from "Edit Line 1" onwards.

```shell 
$ > git reset --hard f8e297f #reset to the point before the revert
$ > git revert --no-commit 0c03424..HEAD
```
This will revert all the commits from "Edit Line 1" onward.  The reason for the ``--no-commit`` parameter in this example, is because _each_ commit will be reverted in turn.  When you use ``--no-commit`` the total reversions are consolidated and you can do a manual commit with an appropriate message. 

### Why is it cool?

Well, undoing disastrous changes is always great.  Both tools, ``git reset`` and ``git revert`` are fantastic in undoing changes.  The majority of the time, you may want to revert a change due to a specific commit, so ``git revert`` is generally preferred, since it only undoes the changes brought about by that specific commit. Furthermore, you can specify the files(s) without reverting changes to other files in the same commit.







Our scenario is a little complicated, but the ``reset`` command has fairly trivial usages too. 



