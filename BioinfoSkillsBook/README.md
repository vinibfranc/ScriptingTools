# Bioinformatics Data Skills

[Book website](http://vincebuffalo.org/book/)

## Chapter 1 -  How to Learn Bioinformatics

1. Biology’s Growing Data
2. New Challenges for Reproducible and Robust Research
    - Reproducibility
    - Sharing data
    - Analyses can lead to findings being susceptible to errors and technical confounding
3. Reproducible Research
4. Robust Research and the Golden Rule of Bioinformatics
    - In wetlab biology, when experiments fail, it can be very apparent, but this is not always true in computing
    - Genomics data also creates its own challenges for robust research. First, most bioinformatics analyses produce intermediate output that is too large and high dimensional to inspect or easily visualize
5. Genomics data also creates its own challenges for robust research. First, most bioinformatics analyses produce intermediate output that is too large and high dimensional to inspect or easily visualize
    - Need to rerun analyses
    - Revisit parts of a project
6. Recommendations for Robust Research
    - Pay Attention to Experimental Design
    - Write Code for Humans, Write Data for Computers
    - Let Your Computer Do the Work For You
    - Make Assertions and Be Loud, in Code and in Your Methods
    - Test Code, or Better Yet, Let Code Test Code
        - How many times is this code called by other code?
        - If this code were wrong how detrimental to the final results would it be?
        - How noticeable would an error be if one occurred?
    - Use Existing Libraries Whenever Possible
    - Treat Data as Read-Only
    - Spend Time Developing Frequently Used Scripts into Tools
    - Let Data Prove That It’s High Quality
7. Recommendations for Reproducible Research
    - Release Your Code and Data
    - Document Everything
    - Make Figures and Statistics the Results of Scripts
    - Use Code as Documentation
    - Continually Improving Your Bioinformatics Data Skills

## Chapter 2 - Setting Up and Managing a Bioinformatics Project

Just as a well-organized laboratory makes a scientist’s life easier, a well-organized and well-documented project makes a bioinformatician’s life easier.

### Project Directories and Directory Structures

Suppose you’re working on SNP calling in maize (Zea mays).

```
$ mkdir zmays-snps
$ cd zmays-snps
$ mkdir data
$ mkdir data/seqs scripts analysis
$ ls -l
```

This is a sensible project layout scheme. Here, data/ contains all raw and intermediate data. As we’ll see, data-processing steps are treated as separate subprojects in this
data/ directory. I keep general project-wide scripts in a scripts/ directory. If scripts contain many files (e.g., multiple Python modules), they should reside in their own subdirectory. Isolating scripts in their own subdirectory also keeps project directories tidy while developing these scripts (when they produce test output files).

Bioinformatics projects contain many smaller analyses—for example, analyzing the quality of your raw sequences, the aligner output, and the final data that will produce figures and tables for a paper.

### Project Documentation

In addition to having a well-organized directory structure, your bioinformatics project also should be well documented. Poor documentation can lead to irreproduci‐ bility and serious errors. The best way to prevent this complexity from causing problems is to document everything extensively.

So what exactly should you document? Here are some ideas:

- Document your methods and workflows: In general, any command that produces results used in your work needs to be documented somewhere.

- Document the origin of all data in your project directory: You need to keep track of where data was downloaded from, who gave it to you, and any other relevant information. “Data” doesn’t just refer to your project’s experimental data—it’s any data that programs use to create output.

- Document when you downloaded data: It’s important to include when the data was downloaded, as the external data source (such as a website or server) might change in the future. For example, a script that downloads data directly from a database might produce different results if rerun after the external database is updated. Consequently, it’s important to document when data came into your repository.

- Record data version information: Many databases have explicit release numbers, version numbers, or names (e.g., TAIR10 version of genome annotation for Arabidopsis thaliana, or Wormbase release WS231 for Caenorhabditis elegans . It’s important to record all version information in your documentation, including minor version numbers.

- Describe how you downloaded the data: For example, did you use MySQL to download a set of genes? Or the UCSC Genome Browser? These details can be useful in tracking down issues like when data is different between collaborators.

- Document the versions of the software that you ran: Good bioinformatics software usually has a command-line option to return the current version. Software managed with a version control system such as Git has explicit identifiers to every version, which can be used to document the precise version you ran. If no version information is available, a release date, link to the software, and download date will suffice.

All of this information is best stored in plain-text README files. Plain text can easily be read, searched, and edited directly from the command line, making it the perfect choice for portable and accessible README files.

```
touch README data/README
```

### Use Directories to Divide Up Your Project into Subprojects

Bioinformatics projects involve many subprojects and subanalyses. For example, the quality of raw experimental data should be assessed and poor quality regions removed before running it through bioinformatics tools like aligners or assemblers. Even before you get to actually analyzing sequences, your project directory can get cluttered with intermediate files.

Creating directories to logically separate subprojects (e.g., sequencing data quality improvement, aligning, analyzing alignment results, etc.) can simplify complex projects and help keep files organized. It also helps reduce the risk of accidentally clobbering a file with a buggy script, as subdirectories help isolate mishaps. Breaking a project down into subprojects and keeping these in separate subdirectories also makes documenting your work easier; each README pertains to the directory it resides in.

### Organizing Data to Automate File Processing Tasks

Because automating file processing tasks is an integral part of bioinformatics, organizing our projects to facilitate this is essential. Organizing data into subdirectories and using clear and consistent file naming schemes is imperative—both of these practices allow us to programmatically refer to files, the first step to automating a task. Doing something programmatically means doing it through code rather than manually, using a method that can effortlessly scale to multiple objects. Programmatically referring to multiple files is easier and safer than typing them all out (because it’s less error prone).

- Shell Expansion Tips

The next command results in the same path organization that we've done before.

```
mkdir -p zmays-snps/{data/seqs,scripts,analysis}
```

- Creating files to see how consistent names help with grammatically working with files

```
$ cd data
$ touch seqs/zmays{A,B,C}_R{1,2}.fastq
$ ls seqs/
zmaysA_R1.fastq zmaysB_R1.fastq zmaysC_R1.fastq
zmaysA_R2.fastq zmaysB_R2.fastq zmaysC_R2.fastq
```

Suppose that we wanted to programmatically retrieve all files that have the sample name zmaysB rather than having to manually specify each file.

```
$ ls seqs/zmaysB*
zmaysB_R1.fastq zmaysB_R2.fastq
```

Suppose a collaborator tells you that the C sample sequences are poor quality, so you’ll have to work with just the A and B samples while C is resequenced. You don’t want to delete zmaysC_R1.fastq and zmaysC_R2.fastq until the new samples
are received, so in the meantime you want to ignore these files. The folks that invented wildcards foresaw problems like this, so they created shell wildcards that allow you to match specific characters or ranges of characters.

```
$ ls zmays[AB]_R1.fastq
zmaysA_R1.fastq zmaysB_R1.fastq
$ ls zmays[A-B]_R1.fastq
zmaysA_R1.fastq zmaysB_R1.fastq
```

- Leading Zeros and Sorting

Another useful trick is to use leading zeros (e.g., file-0021.txt rather than file-21.txt) when naming files. This is useful because lexicographically sorting files (as ls does) leads to the correct ordering.

### Markdown for Project Notebooks

It’s very important to keep a project notebook containing detailed information about the chronology of your computational work, steps you’ve taken, information about why you’ve made decisions, and of course all pertinent information to reproduce your work.

Markdown is just plain-text, which means that it’s portable and programs to edit and read it will exist. Anyone who’s written notes or papers in old versions of word processors is likely familiar with the hassle of trying to share or update out-of-date proprietary formats. For these reasons, Markdown makes for a simple and elegant notebook format.

- Markdown Formatting Basics (see README on chapter2 directory)

```
*emphasis*:  emphasis
**bold**: bold
`inline code`: inline code
<http://website.com/link>: Hyperlink to http://website.com/link
[link text](http://website.com/link): Hyperlink to http://website.com/link, with text “link text”
![alt text](path/to/figure.png): Image, with alternative text
# Header level 1
## Header level 2
### Header level 3
Header level 1
==============
Header level 2
--------------
- Francis Crick
- James D. Watson
- Rosalind Franklin
```

- Using Pandoc to Render Markdown to HTML

Pandoc can convert between a variety of different markup and output formats. Using Pandoc is very simple—to convert from Markdown to HTML, use the --from mark down and --to html options and supply your input file as the last argument:

```
$ pandoc --from markdown --to html notebook.md > output.html
```

## Chapter 3 - Remedial Unix Shell

The Unix shell is the foundational computing environment for bioinformatics. The shell serves as our interface to large bioinformatics programs, as an interactive con‐ sole to inspect data and intermediate results, and as the infrastructure for our pipelines and workflows.

### Why Do We Use Unix in Bioinformatics? Modularity and the Unix Philosophy

We usually don’t think of a bioinformatics project as a “program” but it certainly could be—we could write a single complex program that takes raw data as input, and after hours of data processing, outputs publication figures and a final table of results. For a project like variant calling, this program would include steps for raw sequence read processing, read mapping, variant calling, filtering variant calls, and final data analysis. This program’s code would be massive—easily thousands of lines long.

Unix is the foundational computing environment in bioinformatics because its design philosophy is the antithesis of this inflexible and fragile approach. The Unix shell was designed to allow users to easily build complex programs by interfacing smaller modular programs together.

The Unix shell provides a way for these programs to talk to each other (pipes) and write to and read files (redirection). Unix’s core programs are modular and designed to work well with
other programs. The modular approach of the Unix philosophy has many advantages in bioinformatics:

- With modular workflows, it’s easier to spot errors and figure out where they’re occurring.
- Modular workflows allow us to experiment with alternative methods and approaches, as separate components can be easily swapped out with other components.
- Modular components allow us to choose tools and languages that are appropriate for specific tasks.
- Modular programs are reusable and applicable to many types of data.

The last point to stress about the Unix shell is that it’s incredibly powerful. With simple features like wildcards, it’s trivial to apply a command to hundreds of files, like in the next example:

```
$ rm -rf tmp-data/aligned-reads* # deletes all old large files
$ # versus
$ rm -rf tmp-data/aligned-reads * # deletes your entire current directory
rm: tmp-data/aligned-reads: No such file or directory
```

### Working with Streams and Redirection

Bioinformatics data is often text—for example, the As, Cs, Ts and Gs in sequencing read files or reference genomes, or tab-delimited files of gene coordinates. The text data in bioinformatics is often large, too (gigabytes or more that can’t fit into your computer’s memory at once). This is why Unix’s philosophy of handling text streams is useful in bioinformatics: text streams allow us to do processing on a stream of data rather than holding it all in memory.

1. Redirecting Standard Out to a File

We can combine large files by printing their contents to the standard output stream and redirect this stream from our terminal to the file we wish to save the combined results
to.

Sequence data from: https://github.com/vsbuffalo/bds-files/tree/master/chapter-03-remedial-unix

```
$ cat tb1-protein.fasta
>teosinte-branched-1 protein
LGVPSVKHMFPFCDSSSPMDLPLYQQLQLSPSSPKTDQSSSFYCYPCSPP
FAAADASFPLSYQIGSAAAADATPPQAVINSPDLPVQALMDHAPAPATEL
GACASGAEGSGASLDRAAAAARKDRHSKICTAGGMRDRRMRLSLDVARKF
FALQDMLGFDKASKTVQWLLNTSKSAIQEIMADDASSECVEDGSSSLSVD
GKHNPAEQLGGGGDQKPKGNCRGEGKKPAKASKAAATPKPPRKSANNAHQ
VPDKETRAKARERARERTKEKHRMRWVKLASAIDVEAAAASVPSDRPSSN
NLSHHSSLSMNMPCAAA

$ cat tb1-protein.fasta tga1-protein.fasta
>teosinte-branched-1 protein
LGVPSVKHMFPFCDSSSPMDLPLYQQLQLSPSSPKTDQSSSFYCYPCSPP
FAAADASFPLSYQIGSAAAADATPPQAVINSPDLPVQALMDHAPAPATEL
GACASGAEGSGASLDRAAAAARKDRHSKICTAGGMRDRRMRLSLDVARKF
FALQDMLGFDKASKTVQWLLNTSKSAIQEIMADDASSECVEDGSSSLSVD
GKHNPAEQLGGGGDQKPKGNCRGEGKKPAKASKAAATPKPPRKSANNAHQ
VPDKETRAKARERARERTKEKHRMRWVKLASAIDVEAAAASVPSDRPSSN
NLSHHSSLSMNMPCAAA
>teosinte-glume-architecture-1 protein
DSDCALSLLSAPANSSGIDVSRMVRPTEHVPMAQQPVVPGLQFGSASWFP
RPQASTGGSFVPSCPAAVEGEQQLNAVLGPNDSEVSMNYGGMFHVGGGSG
GGEGSSDGGT
```

To redirect this concatenated sequences to a file, use the operators > or >>.

```
$ cat tb1-protein.fasta tga1-protein.fasta > zea-proteins.fasta
```

2. Redirecting Standard Error

Because many programs use the standard output stream for outputting data, a separate stream is needed for errors, warnings, and messages meant to be read by the user. Standard error is a stream just for this purpose. Like standard output, standard error is by default directed to your terminal. In practice, we often want to redirect the standard error stream to a file so messages, errors, and warnings are logged to a file we can check later.

Trying to access a file not available:

```
$ ls -l tb1-protein.fasta leafy1.fasta
ls: leafy1.fasta: No such file or directory
-rw-r--r-- 1 vinceb staff 0 Feb 21 21:58 tb1.fasta
```

To redirect each stream to separate files, we combine the > operator from the previous section with a new operator for redirecting the standard error stream, 2>:

```
$ ls -l tb1-protein.fasta leafy1.fasta > listing.txt 2> listing.stderr
$ cat listing.txt
-rw-r--r-- 1 vinceb staff 152 Jan 20 21:24 tb1.fasta
$ cat listing.stderr
ls: leafy1.fasta: No such file or directory
```

3. Using Standard Input Redirection

The Unix shell also provides a redirection operator for standard input. Normally standard input comes from your keyboard, but with the < redirection operator you can read standard input directly from a file. Though standard input redirection is less common than >, >>, and 2>, it is still occasionally useful:

```
$ program < inputfile > outputfile
```

In this example, the artificial file inputfile is provided to program through standard input, and all of program’s standard output is redirected to the file outputfile.

### The Almighty Unix Pipe: Speed and Beauty in One

Passing the output of one program directly into the input of another program with pipes is a computationally efficient and simple way to interface Unix programs. This is another reason why bioinformaticians (and software engineers in general) like Unix. Pipes allow us to build larger, more complex tools from smaller modular parts. It doesn’t matter what language a program is written in, either; pipes will work between anything as long as both programs understand the data passed between them.

1. Pipes in Action: Creating Simple Programs with Grep and Pipes

Our pipeline would first remove all header lines (those that begin with >) from the FASTA files, as we only care if sequences have non-nucleotide characters. The remaining sequences of the FASTA file could then be piped to another instance of grep, which would only print lines containing non-nucleotide characters. To make these easier to spot in our terminal, we could also color these matching characters. The entire command would look like:

```
$ grep -v "^>" tb1-protein.fasta | grep --color -i "[^ATCG]"
CCCCAAAGACGGACCAATCCAGCAGCTTCTACTGCTAYCCATGCTCCCCTCCCTTCGCCGCCGCCGACGC
```

2. Combining Pipes and Redirection

Large bioinformatics programs like aligners, assemblers, and SNP callers will often use multiple streams simultaneously. For example, suppose we have two imaginary programs: program1 and program2. Our first program, program1, does some processing on an input file called input.txt and outputs results to the standard output stream and diagnostic messages to the standard error stream. Our second program, program2, takes standard output from program1 as input and processes it. program2 also outputs its own diagnostic messages to standard error, and results to standard output. The tricky part is that we now have two processes outputting to both standard error and standard output. Luckily, we can can combine pipes and redirects easily:

```
$ program1 input.txt 2> program1.stderr | \
  program2 2> program2.stderr > results.txt
```

Occasionally, we need to redirect a standard error stream to standard output. For example, suppose we wanted to use grep to search for “error” in both the standard output and standard error streams of program1.

```
$ program1 2>&1 | grep "error"
```

3. Even More Redirection: A tee in Your Pipe

As mentioned earlier, pipes prevent unnecessary disk writing and reading operations by connecting the standard output of one process to the standard input of another. However, we do occasionally need to write intermediate files to disk in Unix pipelines. These intermediate files can be useful when  debugging a pipeline or when you wish to store intermediate files for steps that take a long time to complete.

```
$ program1 input.txt | tee intermediate-file.txt | program2 > results.txt
```

### Managing and Interacting with Processes

When we run programs through the Unix shell, they become processes until they successfully finish or terminate with an error. There are multiple processes running on your machine simultaneously—for example, system processes, as well as your web browser, email application, bioinformatics programs, and so on. In bioinformatics, we often work with processes that run for a large amount of time, so it’s important we know how to work with and manage processes from the Unix shell.

1. Background Processes

When we type a command in the shell and press Enter, we lose access to that shell prompt for however long the command takes to run. This is fine for short tasks, but waiting for a long-running bioinformatics program to complete before continuing work in the shell would kill our productivity. Rather than running your programs in the foreground (as you do normally when you run a command), the shell also gives you the option to run programs in the background. Running a process in the background frees the prompt so you can continue working.

We can tell the Unix shell to run a program in the background by appending an ampersand (&) to the end of our command. For example:

```
$ program1 input.txt > results.txt &
[1] 26577
```

We also can check what processes we have running in the background with jobs:

```
$ jobs
[1]+ Running program1 input.txt > results.txt
```

To bring a background process into the foreground again, we can use fg (for foreground). fg will bring the most recent process to the foreground. If you have many processes running in the background, they will all appear in the list output by the program jobs. To return a specific background job to the foreground, use fg %num where num is its number in the job list.

```
$ fg
program1 input.txt > results.txt
```

It’s also possible to place a process already running in the foreground into the background. To do this, we first need to suspend the process, and then use the bg command to run it in the background.

```
$ program1 input.txt > results.txt # forgot to append ampersand
$ # enter control-z
[1]+ Stopped program1 input.txt > results.txt
$ bg
[1]+ program1 input.txt > results.txt
```

As with fg earlier, we could also use jobs to see the suspended process’s job ID. If we have multiple running processes, we can specify which one to move to the back‐ ground with bg %num.

There’s a slight gotcha with background processes: although they run in the background and seem disconnected from our terminal, closing our terminal window would cause these processes to be killed.

2. Killing Processes

Occasionally we need to kill a process. It’s not uncommon for a process to demand too many of our computer’s resources or become nonresponsive, requiring that we send a special signal to kill the process. Killing a process ends it for good, and unlike suspending it with a stop signal, it’s unrecoverable. If the process is currently running in your shell, you can kill it by entering Control-C. This only works if this process is running in the foreground, so if it’s in the background you’ll have to use the fg discussed earlier.

3. Exit Status: How to Programmatically Tell Whether Your Command Worked

One concern with long-running processes is that you’re probably not going to wait around to monitor them. How do you know when they complete? How do you know if they successfully finished without an error? Unix programs exit with an exit status,
which indicates whether a program terminated without a problem or with an error. By Unix standards, an exit status of 0 indicates the process ran successfully, and any nonzero status indicates some sort of error has occurred (and hopefully the program prints an understandable error message, too).

The exit status isn’t printed to the terminal, but your shell will set its value to a variable in your shell (aptly named a shell variable) named $?. We can use the echo command to look at this variable’s value after running a command:

```
$ program1 input.txt > results.txt
$ echo $?
0
```

Exit statuses are incredibly useful because they allow us to programmatically chain commands together in the shell. A subsequent command in a chain is run conditionally on the last command’s exit status. The shell provides two operators that implement this: one operator that runs the subsequent command only if the first command completed successfully (&&), and one operator that runs the next command only if the first completed unsuccessfully (||).

It’s best to see these operators in an example. Suppose we wanted to run program1, have it write its output to file, and then have program2 read this output. To avoid the problem of program2 reading a file that’s not complete because program1 terminated with an error, we want to start program2 only after program1 returns a zero (successful) exit code. The shell operator && executes subsequent commands only if previous commands have completed with a nonzero exit status:

```
$ program1 input.txt > intermediate-results.txt && \
program2 intermediate-results.txt > results.txt
```

Using the || operator, we can have the shell execute a command only if the previous command has failed (exited with a nonzero status). This is useful for warning messages:

```
$ program1 input.txt > intermediate-results.txt || \
echo "warning: an error occurred"
```

### Command Substitution

Unix users like to have the Unix shell do work for them — this is why shell expansions like wildcards and brace expansion exist. Another type of useful shell expansion is command substitution. Command substitution runs a Unix command inline and returns the output as a string that can be used in another command. This opens up a lot of useful possibilities.

We use command substitution to run the date program and replace
this command with its output (the string). This is easier to understand through a simpler example:

```
$ grep -c '^>' input.fasta
416
$ echo "There are $(grep -c '^>' input.fasta) entries in my FASTA file."
There are 416 entries in my FASTA file.
```

Using this command substitution approach, we can easily create dated directories using the command date +%F, where the argument +%F simply tells the date program to output the date in a particular format.

```
$ mkdir results-$(date +%F)
$ ls results-2018-04-04
```

Early Unix users were a clever (or lazy) bunch and devised a tool for storing repeated command combinations: alias. Below we have some examples:

```
alias mkpr="mkdir -p {data/seqs,scripts,analysis}"

alias today="date +%F"
mkdir results-$(today)
```

## Chapter 6 - Bioinformatics Data

Data is a requisite of any bioinformatics project. We further our understanding of complex biological systems by refining a large amount of data to a point where we can extract meaning from it. Unfortunately, many tasks that are simple with small or medium-sized datasets are a challenge with the large and complex datasets common in genomics. These challenges include:

- Retrieving data: Whether downloading large sequencing datasets or accessing a web application
hundreds of times to download specific files, retrieving data in bioinformatics can require special tools and skills.

- Ensuring data integrity: Transferring large datasets across networks creates more opportunities for data corruption, which can later lead to incorrect analyses. Consequently, we need to ensure data integrity with tools before continuing with analysis.

- Compression: The data we work with in bioinformatics is large enough that it often needs to be compressed. Consequently, working with compressed data is an essential skill in bioinformatics.

### Retrieving Bioinformatics Data

Suppose you’ve just been told the sequencing for your project has been completed: you have six lanes of Illumina data to download from your sequencing center. Downloading this amount of data through your web browser is not feasible: web browsers are not designed to download such large datasets. Additionally, you’d need to download this sequencing data to your server, not the local workstation where you browse the Internet.

1. Downloading Data with wget and curl

Two common command-line programs for downloading data from the Web are wget and curl.

#### wget

wget is useful for quickly downloading a file from the command line—for example, human chromosome 22 from the GRCh37 (also known as hg19) assembly version:

```
$ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
--2013-06-30 00:15:45-- http://[...]/goldenPath/hg19/chromosomes/chr22.fa.gz
Resolving hgdownload.soe.ucsc.edu... 128.114.119.163
Connecting to hgdownload.soe.ucsc.edu|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 11327826 (11M) [application/x-gzip]
Saving to: ‘chr22.fa.gz’

17% [======> ] 1,989,172 234KB/s eta 66s
```

wget can handle both http and ftp links. One of wget’s strengths is that it can download data recursively. When run with the recursive option (--recursive or -r), wget will also follow and download the pages linked to, and even follow and download links on these pages, and so forth. An exemple is provided below:

```
$ wget --accept "*.gtf" --no-directories --recursive --no-parent http://genomics.someuniversity.edu/labsite/annotation.html
```

#### curl

curl serves a slightly different purpose than wget. wget is great for downloading files via HTTP or FTP and scraping data from a web page using its recursive option. curl behaves similarly, although by default writes the file to standard output. To download
chromosome 22 as we did with wget, we’d use:

```
$ curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz > chr22.fa.gz
```

curl has the advantage that it can transfer files using more protocols than wget, including SFTP (secure FTP) and SCP (secure copy). Curl itself is also a library, meaning in addition to the command-line curl program, Curl’s functionality is wrapped by software libraries like RCurl and pycurl.

### Rsync and Secure Copy (scp)

wget and curl are appropriate for quickly downloading files from the command line, but are not optimal for some heavier-duty tasks. For example, suppose a colleague needs all large sequencing datasets in your project directory that are ignored by Git. A better tool for synchronizing these entire directories across a network is Rsync.

Rsync is a superior option for these types of tasks for a few reasons. First, Rsync is often faster because it only sends the difference between file versions (when a copy already exists or partially exists) and it can compress files during transfers. Second, Rsync has an archive option that preserves links, modification timestamps, permissions, ownership, and other file attributes. This makes Rsync an excellent choice for network backups of entire directories. Rsync also has numerous features and options to handle different backup scenarios, such as what to do if a file exists on the remote host.

rsync’s basic syntax is rsync source destination, where source is the source of the files or directories you’d like to copy, and destination is the destination you’d like to copy these files to. Either source or destination can be a remote host specified in the format user@host:/path/to/directory/.

Let’s look at an example of how we can use rsync to copy over an entire directory to another machine. 

```
$ rsync -avz -e ssh zea_mays/data/ vinceb@[...]:/home/deborah/zea_mays/data
building file list ... done
zmaysA_R1.fastq
zmaysA_R2.fastq
zmaysB_R1.fastq
zmaysB_R2.fastq
zmaysC_R1.fastq
zmaysC_R2.fastq
sent 2861400 bytes received 42 bytes 107978.94 bytes/sec
total size is 8806085 speedup is 3.08
```

Occasionally, we just need to quickly copy a single file over SSH—for tasks where Unix’s cp would be sufficient, but needs to work over an SSH connection. rsync would work, but it’s a bit overkill. Secure copy (scp) is perfect for this purpose. Secure
copy works just like cp, except we need to specify both host and path (using the same user@host:/path/to/file notation as wget). For example, we could transfer a single GTF file to 192.168.237.42:/home/deborah/zea_mays/data/ using:

```
$ scp Zea_mays.AGPv3.20.gtf \ 192.168.237.42:/home/deborah/zea_mays/data/
```

### Data Integrity

Data we download into our project directory is the starting point of all future analyses and conclusions. Although it may seem improbable, the risk of data corruption during transfers is a concern when transferring large datasets. These large files take a long time to transfer, which translates to more opportunities for network connections to drop and bits to be lost. In addition to verifying your transfer finished without error, it’s also important to explicitly check the transferred data’s integrity with checksums. Checksums are very compressed summaries of data, computed in a way that even if just one bit of the data is changed, the checksum will be different.

Data integrity checks are also helpful in keeping track of data versions. In collabora‐ tive projects, our analyses may depend on our colleagues’ intermediate results. When these intermediate results change, all downstream analyses that depend on these results need to be rerun.

### SHA and MD5 Checksums

The two most common checksum algorithms for ensuring data integrity are MD5 and SHA-1. MD5 is an older checksum algorithm, but one that is still commonly used. Both MD5 and SHA-1 behave similarly, but SHA-1 is newer and generally preferred. However, MD5 is more common; it’s likely to be what you encounter if
a server has precomputed checksums on a set of files.

Let’s get acquainted with checksums using SHA-1.

```
$ echo "bioinformatics is fun" | shasum
f9b70d0d1b0a55263f1b012adab6abf572e3030b -
$ echo "bioinformatic is fun" | shasum
e7f33eedcfdc9aef8a9b4fec07e58f0cf292aa67 -
```

We can also use checksums with file input (note that the content of Csyrichta_TAGGACT_L008_R1_001.fastq is fake example data):

```
$ shasum Csyrichta_TAGGACT_L008_R1_001.fastq
fea7d7a582cdfb64915d486ca39da9ebf7ef1d83 Csyrichta_TAGGACT_L008_R1_001.fastq
```

When downloading many files, it can get rather tedious to check each checksum individually. The program shasum has a convenient solution—it can create and validate against a file containing the checksums of files. We can create a SHA-1 checksum file for all FASTQ files in the data/ directory as follows:

```
$ shasum data/*fastq > fastq_checksums.sha
$ cat fastq_checksums.sha
```

Then, we can use shasum’s check option (-c) to validate that these files match the original versions:

```
$ shasum -c fastq_checksums.sha
```

### Looking at Differences Between Data

While checksums are a great method to check if files are different, they don’t tell us how files differ. One approach to this is to compute the diff between two files using the Unix tool diff. Unix’s diff works line by line, and outputs blocks (called hunks) that differ between files.

An example is:

```
$ diff -u gene-1.bed gene-2.bed
```

The option -u tells diff to output in unified diff format, which is a format nearly identical to the one used by git diff.

### Compressing Data and Working with Compressed Data

Data compression, the process of condensing data so that it takes up less space (on disk drives, in memory, or across network transfers), is an indispensable technology in modern bioinformatics. For example, sequences from a recent Illumina HiSeq run when compressed with Gzip take up 21,408,674,240 bytes, which is a bit under 20 gigabytes. Uncompressed, this file is a whopping 63,203,414,514 bytes (around 58 gigabytes). The compression ratio (uncompressed size/compressed size) of this data is approximately 2.95, which translates to a significant
space saving of about 66%.

For the most part, data can remain compressed on the disk throughout processing and analyses. Most well-written bioinformatics tools can work natively with compressed data as input, without requiring us to decompress it to disk first.

#### gzip

The two most common compression systems used on Unix are gzip and bzip2. Both have their advantages: gzip compresses and decompresses data faster than bzip2, but bzip2 has a higher compression ratio (the previously mentioned FASTQ file is only about 16 GB when compressed with bzip2). Generally, gzip is used in bioinformatics to compress most sizable files, while bzip2 is more common for long-term data archiving. We’ll focus primarily on gzip, but bzip2’s tools behave very similarly to gzip.

For example, suppose we have a program that removes low-quality bases from FASTQ files called trimmer (this is an imaginary program). Our trimmer program can handle gzipped input files natively, but writes uncompressed trimmed FASTQ results to standard output. Using gzip, we can compress trimmer’s output in place, before writing to the disk:

```
$ trimmer in.fastq.gz | gzip > out.fastq.gz
```

gzip also can compress files on disk in place.

```
$ ls
in.fastq
$ gzip in.fastq
$ ls
in.fastq.gz
```

Similarly, we can decompress files in place with the command gunzip:

```
$ gunzip in.fastq.gz
$ ls
in.fastq
```

A nice feature of the gzip compression algorithm is that you can concatenate gzip compressed output directly to an existing gzip file, like below.

```
$ ls
in.fastq.gz in2.fastq
$ gzip -c in2.fastq >> in.fastq.gz
```

### Case Study: Reproducibly Downloading Data

- Downloading data:

```
wget ftp://ftp.ensembl.org/pub/release-74/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
```

- Extracting the FASTA headers on this gzipped file:

```
$ zgrep "^>" Mus_musculus.GRCm38.74.dna.toplevel.fa.gz | less
```

- Creating checksum file and comparing them:

```
$ wget ftp://ftp.ensembl.org/pub/release-74/fasta/mus_musculus/dna/CHECKSUMS
$ sum Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
53504 793314
```

- Calculating the SHA-1 sum using shasum:

```
$ shasum Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
01c868e22a981[...]c2154c20ae7899c5f Mus_musculus.GRCm38.74.dna.toplevel.fa.gz
```

- Downloading an accompanying GTF from Ensembl and the CHECKSUMS file for this directory:

```
$ wget ftp://ftp.ensembl.org/pub/release-74/gtf/mus_musculus/Mus_musculus.GRCm38.74.gtf.gz
$ wget ftp://ftp.ensembl.org/pub/release-74/gtf/mus_musculus/CHECKSUMS
```

- Ensuring that our checksums match those in the CHECKSUMS file:

```
$ sum Mus_musculus.GRCm38.74.gtf.gz
00985 15074
$ shasum cf5bb5f8bda2803410bb04b708bff59cb575e379 Mus_musculus.GRCm38.74.gtf.gz
```

## Chapter 7
## Chapter 8
## Chapter 9
## Chapter 10
## Chapter 11
## Chapter 12
## Chapter 13