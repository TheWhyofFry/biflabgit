#Ignore this guide in its current form
# Git guide 2 - Resolving conflicts

In an ideal world, everyone involved in a group porject will edited only that which was delegated to them, with the project being modular (split into enough bits) so that there are no more than one person working on a project. In reality, this is a little more tricky and the inevitable conflict will result.  

You may have thought about what we've dealt with previously and reckon: 

> "All well and good, sir, but what if my friend and I work on the same file? Will ``Git`` know who did what and change things in a rational way?"

Well.  Partly yes, but mostly . . no. To give a simple example, lets consider the following block of text - some lyrics from the Pink Floyd song, "Eclipse".  Imagine that only a partial copy of the lyrics are available, because ... reasons.

The copy you receive cannot possibly be complete, because there are only five lines included! How rude. You could listen to the song and capture the rest of the lyrics manually.  Instead, you take this opportunity to test out some Git skills you've learned over the last couple of weeks.  You and a friend (Bob) contribute to correcting and completing the lyrics.  You create an initial commit with the first version of the lyrics (Commit1). Both you and your friend pull the same version of the repository and both of you are in sync.  The initial commit, ``Commit1``.  Bob and you edit your respective copies of ``commit1.txt`` and commit the changes. Bob ``pushes`` his changes to the repository before you, so when you try to do it, ``Git`` complains about a conflict.  

To resolve these conflicts, it helps to use a tool that can show the changes between _three_ copies of the file.  You may think that you only need to look at differences between your file and Bob's file, but remember that you each made changes to the same file. You need to see where the changes are:

 * Inconsequential - you both edited separate parts of the file, there are no lines edited by either you or Bob that are 
 * Conflicting - Both you and Bob made changes to a line(s) of the lyrics that are _not_ the same and needs to be resolved

To illustrate this, we'll use three copies of imaginary commits (in the data/guide2 folder).  

 * ``commit1.txt`` - The first commit that both you and Bob received
 * ``bob.txt`` - Changes made to ``commit1.txt`` by *Bob*.
 * ``you.txt`` - Changes made to ``commit1.txt`` by *you*.


First, download the tool "KDiff3" [here](https://sourceforge.net/projects/kdiff3/files/latest/download). This tool will help with the visualization and resolving conflicts in the edit.  Download the respective files, ``commit1.txt``, ``bob.txt`` and ``you.txt`` to a location of your choosing. Open up a terminal and navigate to the folder containing these files - preferably have a new folder that does not contain any other files.  You can also clone this guide's repository and navigate to the ``data/guide2`` folder.  From the command line (if you're using bash):

```shell
$ biflabgit/data/guide2> kdiff3 commit1.txt bob.txt you.txt
```

The first parameter is the ``base`` file, the copy of the file right at the beginning that you and Bob both had. The second is Bob's edit and the third is yours.  The order is important, since it is _your_ file that will contain the conflict resolved version in the end. After ``KDiff3`` starts up it will present you with three columns. 

![kdiff3](images/guide3_kdiff3.png)

The column order follows the order of the file parameters you've given.  The job here is for you to _merge_ the changes into one file, a final version of the edit.  In the majority of cases, the merges can be auto-resolved.  Click on the right-most column (the one for ``you.txt``). Click on the ``Merge -> Merge Current file`` menu option.  ``KDiff3`` will report that there are four conflicts in the file. It could auto-resolve three out of four. 

Instead of doing things line-by-line, diff tries to get blocks of text that are different between files.  The first block deals with the lines:

> everyone you meat
> All you dislike

These lines have both been changed by you and Bob to:

> everyone you *meet*
> All *that* you *slight*

Since both you and Bob's edits are the same, there is no conflict and this part can be auto-merged. The next block of changes are insertions.  These lines did not exist in the original ``commit1.txt`` file and are also shared between your edit and Bob's edit.  So Bob added the lines (for this block):

> All that is now
> All that is gone
> All that's to come

The same block is hilighted in your edit.  However, there is a slight difference between the first line of the block.

> All *that is* now (Bob)
> All *that's* now (You)

This conflict cannot be auto-resolved. When you initially chose the ``Merge -> Merge Current File`` option in ``KDiff3``, it automatically moves to the first line where it cannot automatically resolve a conflict. In the top-right of the window, there are three buttons "A", "B" and "C". Each of them represents the files of column A, B and C, read left to right. If you decide Bob's version is correct, you will click on B.  If you think your version is correct, you will click on C.  Note, however, that you can both select and deselect multiple columns.  However, to resolve the merge (where there is unresolvable conflict) you need to choose _only_ the column of the change you accept.  Bob's version is correct, so choose "B". 

The last block added by you wasn't added by Bob, so the change can be left at the automatic of "C".  You can save the changes, by chosing the ``File -> Save`` option.  It will save the changes to *``you.txt``* unless you specify otherwise.  By default, a copy of the original ``you.txt`` will be made as ``you.txt.orig`` to save you some heartache ;-).  


## Diff 

Another way to see what's going on behind the accepted changes to you, let's see how the new copy of ``you.txt`` (the output from ``KDiff3``) compares to ``bob.txt`` and then how ``bob.txt`` compares with ``you.txt.orig``.  The tool to use is ``diff``.  Diff works similarly to ``KDiff3``, but is focused on comparing two files (three files can be done with ``diff3``). Diff takes in two filenames, compares them and outputs a patch string.  This patch string can be applied to the first file to produce the second file.  Behind the scenes, this is how ``Git`` compares variations of files between commits. 

Using the ``diff`` tool:
```shell
$ biflabgit/data/guide2 > diff bob.txt you.txt
```
```diff
23a24,27
> And everything under the sun is in tune
> But the sun is eclipsed by the moon.
> There is no dark side of the moon really.
> Matter of fact it's all dark.
```

Compared to ``bob.txt`` the file ``you.txt`` adds four lines from line 24 to 27, indicated by ``a``. If you compare ``bob.txt`` to ``you.txt.orig``, diff will output the following:

```diff
21c21
< All that is now
---
> All that's now
23a24,27
> And everything under the sun is in tune
> But the sun is eclipsed by the moon.
> There is no dark side of the moon really.
> Matter of fact it's all dark.
```

So the diff output still shows the added lines, but it also shows that there is a change in line 21. When we merged the files with ``KDiff3`` the conflict in line 21 was resolved, so the only difference between ``bob.txt`` and ``you.txt`` were the added lines.  



