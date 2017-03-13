# About this guide

This guide is a very brief overview of ``git`` and will be alinged with the lectures as presented in the lab.  It is by no means a complete guide. To recap, ``git`` is a version control system. This means that instead of holding single versions of a file, or having copies of files with horrible names such as ``final_report_i_mean_it_this_time.tex`` a history is held as the original file ``report.tex`` changes over time.  By keeping track of _all_ the changes made not only to specific files, but to the entire directory (folder) of your work, project management becomes simpler. Before entering the *_gittening_* we'll revisit a traditional folder setup. Consider the following folder.

 
```shell
werner@linux-dcgb:/tmp/test> ls -l
total 316
-rw-r--r-- 1 werner users 300317 Mar 13 20:49 DATA.txt
-rw-r--r-- 1 werner users  15149 Mar 13 20:49 README.txt
-rw-r--r-- 1 werner users   3125 Mar 13 20:49 REPORT.txt
werner@linux-dcgb:/tmp/test> 
```

There is implicitly no information on previous versions of this file. Your supervisor is reading the introduction of your ``REPORT.TXT`` and emails you a copy back with corrections, ``REPORT_SS.TXT``. This is now the edited and annotated version of the ``REPORT.TXT``. Assuming your supervisor is correct, you rename ``REPORT.TXT`` to ``REPORT_DRAFT1.TXT`` and ``REPORT_SS.TXT`` to ``REPORT_LATEST.TXT``. It's easy to imagine how confusing this process becomes at the 10th iteration of this process. 


```shell
werner@linux-dcgb:/tmp/test> ls -l
-rw-r--r-- 1 werner users  15149 Mar 13 20:49 README_SS.txt
-rw-r--r-- 1 werner users  15149 Mar 29 20:49 README_LATEST.txt
-rw-r--r-- 1 werner users  15149 Mar 15 20:49 README_DRAFT1.txt
-rw-r--r-- 1 werner users  15149 Mar 13 20:49 README.txt
-rw-r--r-- 1 werner users  15149 Mar 13 20:49 README_SS.txt
-rw-r--r-- 1 werner users  15149 Mar 29 20:49 README_LATEST.bak
-rw-r--r-- 1 werner users  15149 Mar 15 20:49 README_DRAFT1.txt
-rw-r--r-- 1 werner users  15149 Mar 13 20:49 README.txt
werner@linux-dcgb:/tmp/test> 
```



Wouldn't it be great if we had intermediate copies of files where changes can be implicitly reviewed?  You tentatively introduce your supervisor to Git. You get the task of re-organizing the fruits of your research (bad code and all). Enthusiastically, you look up git and decide to get straight to work. Unfortunately, you quickly realize that even though this tool is purpoted to be really really _*really*_ awesome, its use (initially) is not trivial. Add? Push? Pull? Commit? Clone? Blame?! You correctly conclude that it would be better to actually learn git from scratch using toy folders as a starting point. 


## Almost there

A few more things before we start. This guide is rendered in a neat, pretty way as you're reading it.  It is written in "markdown".  You may be familiar with "markup" languages such as HTTML for rendering pretty pages in a browser.  Markdown is a lightweight alternative that allows for the creation of professional looking documents without too much hassle.  Oh, and because the strucutre is simple, it is a trivial task to look at ``diffs`` (more on this later). 

You're viewing this document on GitHub. GitHub is an online resource to store your repositories. The advantage of this is obvious, as it becomes easier toe collaborate with users from all over the world to your project. 



## Creating a git repository

It isn't absolutely crucial to create a repository on GitHub, it is prudent to understand the workings of git in the context of a remote repository as you may very well be required to work on projects hosted on Github in the (possibly near :-) ) future.

So first, create an account at [GitHub](https://github.com). The steps are rather intuitive.  Once you've created your account - and please take your time, there is no rush - create an empty repository called "simplerepo". 


![Repo add](images/guide1_newrepo.png)
![Repo add](images/guide1_name.png)

![Repo creation](images/guide1_emptyrepo.png)

##Fetching the repository from GitHub
So let's fetch this repo. We can open a terminal in your favourite OS.  Invoke the following command.

```shell
$:~/Sources/git> git pull https://github.com/[yourchosenusername]/simplerepo.
Cloning into 'simplerepo'...
warning: You appear to have cloned an empty repository.
Checking connectivity... done.
```

Since there is nothing in the repository, ``git`` will let you know that there is nothing in the repository yet.  So let's create a file.  You might want to let people know what this particular repository is all about.  By default, GitHub will display the contents of a ``README.md`` file.  The extension ".md" indicates that the file is in "markdown" format - like this guide. For now, we'll just include a simple line of text to ``README.md``.

```shell
$:~/Sources/git> cd simplerepo

$:~/Sources/git/simplerepo> echo "#This is the README of simplerepo." >> README.md

$:~/Sources/git/simplerepo> cat README.md

#This is the README of simplerepo.
```


##Adding files to the staging area
We need to add this file to the repository. The ``git add`` command adds a file to the _staging_ area.  You can think of this as a state of "limbo".  We can review this ``add`` with the ``git status`` command.

```shell
$:~/Sources/git/simplerepo> git add README.md 
$:~/Sources/git/simplerepo> git status
On branch master

Initial commit

Changes to be committed:
  (use "git rm --cached <file>..." to unstage)

        new file:   README.md
```

There are several components to this message.


* ``On branch master`` : We can have multiple versions of our repository. Generally, we isolate changes to a different _branch_ (Think of a species split) than the master (e.g. a development branch) and when we want to integrate these changes, we _merge_ them into the main branch.
* *``Initial commit``* : We haven't comitted anything to the repository yet
* *``new file:   README.md``* : The file we just added to the staging area


##Commiting the changes to the repository
In order to fully integrate this file with the repository (making the changes for realsies) we need to issue the ``git commit``. The ``commit`` will _commit_ these changes we've made to the _repository_. 

```shell
$:~/Sources/git/simplerepo> git commit -m "Added a simple line to README.md"
[master (root-commit) c464535] Added a simple line to README.md
 1 file changed, 1 insertion(+)
 create mode 100644 README.md
```

If you go to the repository on GitHub, you'll notice that it is still empty.  We've only comitted the changes we've made to the _local_ copy of our repository.  To make changes to the _remote_ repository, the one on GitHub, ``git push`` is used. It is also necessary to add the _origin_ of this repository so it knows where to ``push`` the commit to. However, because we cloned the repository from GitHub, the remote _origin_ is already there.  For the purpose of this exercise we remove our _origin_ from the repository.

```
$:~/Sources/git/simplerepo> git remote rm origin
```

##Adding a remote repository

The ``git remote add`` command is used to add a remote repository where the changes will be _pushed_ to. 

```
$:~/Sources/git/simplerepo> git remote add origin https://github.com/[yourusername]/simplerepo.git
$:~/Sources/git/simplerepo> git remote 
origin
```

##Pushing the changes to the remote repository 
Great, now we've added a remote repository. The changes can be pushed to the remote github.  The ``git push`` command needs a _origin_ and _branch_ of what needs to be committed. It is also possible to invoke the ``-u`` parameter after ``push`` to save the default _remote orgin_ so that with future commits, just using ``git push`` will suffice.

```
$:~/Sources/git/simplerepo> git push -u origin master
Counting objects: 3, done.
Writing objects: 100% (3/3), 299 bytes | 0 bytes/s, done.
Total 3 (delta 0), reused 0 (delta 0)
To https://github.com/thewhyoffry/simplerepo.git
   51ac36b..53a536e  master -> master
Branch master set up to track remote branch master from origin.
```

The funny numbers are the leading values of a _hash_ calculated for this commit.  Simply put, hashses are unique IDs so we have reference points when we want to refer to a specific commit. The value is actually 40 characters long, but for simplicity only the first seven characters of the has is given.  

![Remote initial commit](images/guide1_initialcommit.png)


