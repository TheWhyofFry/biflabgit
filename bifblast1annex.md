# Tips to get started with BioPython


## Pakcage installs on your normal PCs


First, you need to add an environment variable for your BASH (Shell). Edit the file ``~/.bashrc`` and add the following lines to it:

```bash
export proxystring=username:password@intacc.up.ac.za:8080
export http_proxy=http://$proxystring
export https_proxy=https://$proxystring
```

Where ``username`` and ``password`` are your internet username and passwords. These settings will only take effect when you open a new terminal.  If you wish for the settings to work in the current terminal too, execute the commands in that terminal. 

I understand that most of you do not know the root passwords to your machines.  I have told H. the password.  Do NOT write this password down on a piece of paper. 

I don't use CentOS, but as far as I know, you can use ``yum`` to install packages (central software packages).  So, to install ``git``, run:

```shell
$ > sudo -E yum install git
```

The command ``sudo`` will temporarily elevate your user's privileges and allow you to install system packages. 


Similarly, with the added proxy settings, you can now install python packages.