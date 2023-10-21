---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Setting up the environment

## 1. Folder structure

First, we will have a look at the folder structure we will use during the workshop. For any project, it is important to organise your various files into subfolders. Through the course you will be navigating between folders to run your analyses. It may seem confusing or tedious at first, especially as there are only a few files we will be dealing with during this workshop. However, remember that for many projects you can easily generate hundreds of files, so it is best to start with good practices from the beginning.

We will be using the following folder structure for this experiment:

```none
└── sequenceData
    ├── 0-metadata   
    ├── 1-scripts      
    ├── 2-raw    
    ├── 3-fastqc
    ├── 4-demux
    ├── 5-filter
    ├── 6-quality
    ├── 7-refdb
    └── 8-final
```

To set up this folder structure, let us navigate to the starting directory of this workshop. We can do this using the `cd` command, which stands for *change directory*.

```{warning}
The code below assumes you are in the home folder on the HKU supercomputer. Make sure to alter the path that specifies the location of the starting folder on your system if this is not your starting point.
```

```{code-block} bash
cd ednaw01/
```

It is always a good idea to check if the code executed as expected. In this case, we can verify we are in the correct working directory by using the `pwd` command, which stands for *print working directory*.

```{code-block} bash
pwd
```

````{admonition} Output
:class: note, dropdown
```
/home/ednaw01/ednaw01
```
````

We can also list all the documents and subfolders in our working directory using the `ls` command. Additionally, we will specify the `-ltr` parameters to the `ls` command. The `l` parameter is providing a list with long format that includes the permissions, which will come in handy when we are creating our scripts (more on that later). The `t` parameter is sorting the output by time and date, while the `r` parameter is printing the list in reverse order, i.e, the newest files will be printed as the latest in the list.

```{code-block} bash
ls -ltr
```

````{admonition} Output
:class: note, dropdown
```
-rw-r--r-- 1 ednaw01 others     12097 Oct 18 10:31 metadata-COI-selected-updated.txt
-rw-r--r-- 1 ednaw01 others 112435863 Oct 18 10:36 newDBcrabsCOIsintax.fasta
-rw-r--r-- 1 ednaw01 others 568074346 Oct 18 10:40 HKeDNAworkshop2023.zip
-rw-r--r-- 1 ednaw01 others       169 Oct 18 11:36 module_load
-rw-r--r-- 1 ednaw01 others       776 Oct 18 11:37 example-script.sh
```
````

The output of the `ls` command provides us with a list of the documents in our starting folder, which are the documents we need to complete the tutorial. More on this in the next section.

Before we go more in detail about the starting files, let us first set up the abovementioned folder structure. Within the command line interface (CLI), we can use the `mkdir` command to set up all the folders and subfolders in one go. Additionally, we will use the `-p` parameter to automatically make any necessary parent directories that might not yet exist.

```{code-block} bash
mkdir -p sequenceData/0-metadata sequenceData/1-scripts sequenceData/2-raw sequenceData/3-fastqc sequenceData/4-demux sequenceData/5-filter sequenceData/6-quality sequenceData/7-refdb sequenceData/8-final
```

When we now use the `ls -ltr` command again, we see that the directory **sequenceData** has been added to the output list.

```{code-block} bash
ls -ltr
```

````{admonition} Output
:class: note, dropdown
```
-rw-r--r--  1 ednaw01 others       776 Oct 21 10:52 example-script.sh
-rw-r--r--  1 ednaw01 others 568074346 Oct 21 10:53 HKeDNAworkshop2023.zip
-rw-r--r--  1 ednaw01 others     12097 Oct 21 10:53 metadata-COI-selected-updated.txt
-rw-r--r--  1 ednaw01 others       169 Oct 21 10:53 module_load
-rw-r--r--  1 ednaw01 others 112435863 Oct 21 10:53 newDBcrabsCOIsintax.fasta
drwxr-xr-x 11 ednaw01 others        11 Oct 21 10:53 sequenceData
```
````

Note that the `sequenceData` directory has appeared last in the printed list due to the `r` parameter we passed to the `ls` command, as it was the most recently created. When we'll be going through our bioinformatic pipeline and create multiple files at various stages, the `r` parameter will make it easy to see which files were recently created, as they will appear at the bottom and will avoid us having to scroll through the list of files.

Also note the `d` at the beginning of the line of `sequenceData`. These are the permissions that the `l` parameter that we passed to the `ls` command has given us access to. Several other letters are shown for `sequenceData` and other files as well. Below is a list of what they mean:

1. `d`: directory - indicate that this is in fact a folder or directory, rather than a file.
2. `r`: read - specifies that we have access to read the file.
3. `w`: write - specifies that we can write to the file.
4. `x`: execute - specifies that we have permission to execute the file.
5. `@`: extended - a novel symbol for MacOS indicating that the file has extended attributes (MacOS specific).

## 2. Starting files

For the bioinformatic and statistical analysis, we need several starting files, including a zipped file containing all the sequencing data (**HKeDNAworkshop2023.zip**), a sample metadata file for the statistical analysis (**metadata-COI-selected-updated.txt**), our reference database for the taxonomy assignment (**newDBcrabsCOIintax.fasta**), a template script file where we will paste our code to run on the supercomputer (**example-script.sh**), and a script file that automatically loads all the necessary software on the supercomputer. The last two files are only necessary when working on the HKU supercomputer.

```{admonition} Starting files
:class: important
1. HKeDNAworkshop2023.zip
2. metadata-COI-selected-updated.txt
3. newDBcrabsCOIsintax.fasta
4. example-script.sh
5. module_load
```

With the folder structure set up, let's move the starting files to their respective subfolders using the `cp` command.

```{code-block} bash
cp HKeDNAworkshop2023.zip sequenceData/2-raw
cp metadata-COI-selected-updated.txt sequenceData/0-metadata
cp newDBcrabsCOIsintax.fasta sequenceData/7-refdb
cp example-script.sh sequenceData/1-scripts
cp module_load sequenceData/1-scripts
```

Again, we can check if the command executed as expected by listing the files within the subfolder **2-raw** using the `ls -ltr` command. Note that we do not need to first use the `cd` command to move to the directory for which we would like to list all of the files, but that we can specify which directory we would like to list by referring to it after the `ls -ltr` command.

```{code-block} bash
ls -ltr sequenceData/2-raw
```

````{admonition} Output
:class: note, dropdown
```
-rw-r--r-- 1 ednaw01 others 568074346 Oct 21 10:59 HKeDNAworkshop2023.zip
```
````

````{admonition} Exercise 1
:class: hint

Moving files around using the Terminal can also be accomplished through the `mv` command. What is the difference between `mv` and `cp`? Why is it better to use `cp`?

```{admonition} Answer 1
:class: title, dropdown
The `mv` command will move the files from one folder to another. Hence, the files will be removed from the initial directory. If something goes wrong during the process or during the bioinformatic pipeline, such actions cannot be reverted back and data will be lost. 

With the current setup, if anything will go wrong during our pipeline, we can simply remove everything within the **sequenceData** directory using the `rm -r sequenceData*` command and start over, since we still have all the information and files within our initial starting folder.
```
````

**Now that we have everything set up, we can get started with the bioinformatic processing of our data!**
