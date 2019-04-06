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

- **Document your methods and workflows**: In general, any command that produces results used in your work needs to be documented somewhere.

- **Document the origin of all data in your project directory**: You need to keep track of where data was downloaded from, who gave it to you, and any other relevant information. “Data” doesn’t just refer to your project’s experimental data—it’s any data that programs use to create output.

- **Document when you downloaded data**: It’s important to include when the data was downloaded, as the external data source (such as a website or server) might change in the future. For example, a script that downloads data directly from a database might produce different results if rerun after the external database is updated. Consequently, it’s important to document when data came into your repository.

- **Record data version information**: Many databases have explicit release numbers, version numbers, or names (e.g., TAIR10 version of genome annotation for Arabidopsis thaliana, or Wormbase release WS231 for Caenorhabditis elegans . It’s important to record all version information in your documentation, including minor version numbers.

- **Describe how you downloaded the data**: For example, did you use MySQL to download a set of genes? Or the UCSC Genome Browser? These details can be useful in tracking down issues like when data is different between collaborators.

- **Document the versions of the software that you ran**: Good bioinformatics software usually has a command-line option to return the current version. Software managed with a version control system such as Git has explicit identifiers to every version, which can be used to document the precise version you ran. If no version information is available, a release date, link to the software, and download date will suffice.

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

## Chapter 7 - Unix Data Tools

In this chapter, we’ll see how we can combine the Unix shell with command-line data tools to explore and manipulate data quickly. Understanding how to use Unix data tools in bioinformatics isn’t only about learning what each tool does, it’s about mastering the practice of connecting tools together — creating programs from Unix pipelines. By connecting data tools together with pipes, we can construct programs that parse, manipulate, and summarize data.

### When to Use the Unix Pipeline Approach and How to Use It Safely

Many tasks in bioinformatics are of this nature: we want to get a quick answer and keep moving forward with our project. We could write a custom script, but for simple tasks this might be overkill and would take more time than necessary.

For larger, more complex tasks it’s often preferable to write a custom script in a language like Python (or R if the work involves lots of data analysis). While shell approaches (whether a one-liner or a shell script) are useful, these don’t allow for the
same level of flexibility in checking input data, structuring programs, use of data structures, code documentation, and adding assert statements and tests as languages like Python and R.

### Inspecting and Manipulating Text Data with Unix Tools

Tabular plain-text data formats are used extensively in computing. The basic format is incredibly simple: each row (also known as a record) is kept on its own line, and each column (also known as a field) is separated by some delimiter. There are three flavors you will encounter: tab-delimited, comma-separated and variable space-delimited.

In this chapter, we’ll work with very simple genomic feature formats: BED (three-column) and GTF files. These file formats store the positions of features such as genes, exons, and variants in tab-delimited format.

1. Inspecting Data with Head and Tail

Taking a look at the top of a file:

```
$ head Mus_musculus.GRCm38.75_chr1.bed
1 3054233 3054733
1 3054233 3054733
1 3054233 3054733
1 3102016 3102125
1 3102016 3102125
1 3102016 3102125
1 3205901 3671498
1 3205901 3216344
```

We can also control how many lines we see with head through the -n argument:

```
$ tail -n 3 Mus_musculus.GRCm38.75_chr1.bed
1 195240910 195241007
1 195240910 195241007
1 195240910 195241007
```

Sometimes it’s useful to see both the beginning and end of a file—for example, if we have a sorted BED file and we want to see the positions of the first feature and last feature.

```
$ (head -n 2; tail -n 2) < Mus_musculus.GRCm38.75_chr1.bed
1 3054233 3054733
1 3054233 3054733
1 195240910 195241007
1 195240910 195241007
```

We can also search for strings and get the first occurence with head:

```
$ grep 'gene_id "ENSMUSG00000025907"' Mus_musculus.GRCm38.75_chr1.gtf | head -n 1
1 protein_coding gene 6206197 6276648 [...] gene_id "ENSMUSG00000025907" [...]
```

Under the hood, your shell sends a signal to other programs in the pipe called SIGPIPE—much like the signal that’s sent when you press Control-c (that signal is SIGINT). When building complex pipelines
that process large amounts of data, this is extremely important. It means that in a pipeline like:

```
$ grep "some_string" huge_file.txt | program1 | program2 | head -n 5
```

#### less

less is also a useful program for a inspecting files and the output of pipes. less is a terminal pager, a program that allows us to view large amounts of text in our terminals.

less runs more like an application than a command: once we start less, it will stay open until we quit it. Let’s review an example—in this chapter’s directory in the book’s GitHub repository, there’s a file called contaminated.fastq. Let’s look at this with less:

```
$ less contaminated.fastq
```

less is also crucial when iteratively building up a pipeline—which is the best way to construct pipelines. Suppose we have an imaginary pipeline that involves three programs, step1, step2, and step3. Our finished pipeline will look like ```step1 input.txt | step2 | step3 > output.txt```. However, we want to build this up in pieces, running step1 input.txt first and checking its output, then adding in step3 and checking that output, and so forth. The natural way to do this is with less:

```
$ step1 input.txt | less # inspect output in less
$ step1 input.txt | step2 | less
$ step1 input.txt | step2 | step3 | less
```

### Plain-Text Data Summary Information with wc, ls, and awk

In addition to peeking at a file with head, tail, or less, we may want other bits of summary information about a plain-text data file like the number of rows or columns. With plain-text data formats like tab-delimited and CSV files, the number of rows is usually the number of lines. We can retrieve this with the program wc (for word count):

```
$ wc Mus_musculus.GRCm38.75_chr1.bed
81226 243678 1698545 Mus_musculus.GRCm38.75_chr1.bed
```

By default, wc outputs the number of words, lines, and characters of the supplied file. It can also work with many files:

```
$ wc Mus_musculus.GRCm38.75_chr1.bed Mus_musculus.GRCm38.75_chr1.gtf
81226 243678 1698545 Mus_musculus.GRCm38.75_chr1.bed
81231 2385570 26607149 Mus_musculus.GRCm38.75_chr1.gtf
162457 2629248 28305694 total
```

Often, we only care about the number of lines. We can use option -l to just return the number of lines:

```
$ wc -l Mus_musculus.GRCm38.75_chr1.bed
81226 Mus_musculus.GRCm38.75_chr1.bed
```

You might have noticed a discrepancy between the BED file and the GTF file for this chromosome 1 mouse annotation. To inspect this, we can use head:

```
$ head -n 5 Mus_musculus.GRCm38.75_chr1.gtf
```

Another bit of information we usually want about a file is its size. The easiest way to do this is with our old Unix friend, ls, with the -l option:

```
$ ls -l Mus_musculus.GRCm38.75_chr1.bed
-rw-r--r-- 1 vinceb staff 1698545 Jul 14 22:40 Mus_musculus.GRCm38.75_chr1.bed
```

There’s one other bit of information we often want about a file: how many columns it contains. Let’s use an awk one-liner to return how many fields a file contains:

```
$ awk -F "\t" '{print NF; exit}' Mus_musculus.GRCm38.75_chr1.bed
3
```

awk was designed for tabular plain-text data processing, and consequently has a built-in variable NF set to the number of fields of the current dataset. This simple awk one-liner simply prints the number of fields of the first row of the Mus_musculus.GRCm38.75_chr1.bed file, and then exits.

To see how many columns of data there are, we need to first chop off the comments and then pass the results to our awk one-liner. One way to do this is with a tail trick we saw earlier:

```
$ tail -n +5 Mus_musculus.GRCm38.75_chr1.gtf | head -n 1
#!genebuild-last-updated 2013-09
$ tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | head
1 pseudogene gene 3054233 3054733 . + . [...]
$ tail -n +6 Mus_musculus.GRCm38.75_chr1.gtf | awk -F "\t" '{print NF; exit}'
16
```

Otherwise, using the program grep (which we’ll talk more about in, we can easily exclude lines that begin with “#”:

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | head -n 3
1 pseudogene gene 3054233 3054733 . + . [...]
1 unprocessed_pseudogene transcript 3054233 3054733 . + . [...]
1 unprocessed_pseudogene exon 3054233 3054733 . + . [...]
```

### Working with Column Data with cut and Columns

When working with plain-text tabular data formats like tab-delimited and CSV files, we often need to extract specific columns from the original file or stream. For example, suppose we wanted to extract only the start positions (the second column) of the
Mus_musculus.GRCm38.75_chr1.bed file. The simplest way to do this is with cut.

```
$ cut -f 2 Mus_musculus.GRCm38.75_chr1.bed | head -n 3
3054233
3054233
3054233
```

Using cut, we can convert our GTF for Mus_musculus.GRCm38.75_chr1.gtf to a three-column tab-delimited file of genomic ranges (e.g., chromosome, start, and end position). We’ll chop off the metadata rows using the grep command covered earlier, and then use cut to extract the first, fourth, and fifth columns (chromosome, start, end):

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1,4,5 | head -n 3
1 3054233 3054733
1 3054233 3054733
1 3054233 3054733
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1,4,5 > test.txt
```

cut also allows us to specify the column delimiter character. So, if we were to come across a CSV file containing chromosome names, start positions, and end positions, we could select columns from it, too:

```
$ head -n 3 Mus_musculus.GRCm38.75_chr1_bed.csv
1,3054233,3054733
1,3054233,3054733
1,3054233,3054733
$ cut -d, -f2,3 Mus_musculus.GRCm38.75_chr1_bed.csv | head -n 3
3054233,3054733
3054233,3054733
3054233,3054733
```

### Formatting Tabular Data with column

As you may have noticed when working with tab-delimited files, it’s not always easy to see which elements belong to a particular column. For example:

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f1-8 | head -n3
1 pseudogene gene 3054233 3054733 . + .
1 unprocessed_pseudogene transcript 3054233 3054733 . + .
1 unprocessed_pseudogene exon 3054233 3054733 . + .
```

column -t produces neat columns that are much easier to read:

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f 1-8 | column -t
| head -n 3
1 pseudogene gene 3054233 3054733 . + .
1 unprocessed_pseudogene transcript 3054233 3054733 . + .
1 unprocessed_pseudogene exon 3054233 3054733 . + .
```

In general, it’s easier to make computer-readable data attractive to humans than it is to make data in a human-friendly format readable to a computer. Unfortunately, data in formats that prioritize human readability over computer readability still linger in
bioinformatics.

### The All-Powerful Grep

grep is one of the most powerful Unix data tools. First, it’s important to mention grep is fast. Really fast. If you need to find a pattern (fixed string or regular expression) in a file, grep will be faster than anything you could write in Python. This demonstrates a point: if computational speed is our foremost priority (and there are many cases when it isn’t as important as we think), Unix tools tuned to do certain tasks really often are the fastest implementation.

grep requires two arguments: the pattern (the string or basic regular expression you want to search for), and the file (or files) to search for it in. For example, if we want to find a gene in a file, we can do the following:

```
$ grep "Olfr418-ps1" Mus_musculus.GRCm38.75_chr1_genes.txt
ENSMUSG00000049605 Olfr418-ps1
```

One useful option when using grep is --color=auto. This option enables terminal colors, so the matching part of the pattern is colored in your terminal.

For example, suppose you wanted a list of all
genes that contain “Olfr,” except “Olfr1413.” Using -v and chaining together to calls to
grep with pipes, we could use:

```
$ grep Olfr Mus_musculus.GRCm38.75_chr1_genes.txt | grep -v Olfr1413
```

But this command would also exclude genes like
“Olfr1413a” and “Olfr14130.” “Olfr14130.” But we can get around this by using -w, which matches entire words (surrounded by whitespace). Let’s look at how this works with a simpler toy example:

```
$ cat example.txt
bio
bioinfo
bioinformatics
computational biology
$ grep -v bioinfo example.txt
bio
computational biology
$ grep -v -w bioinfo example.txt
bio
bioinformatics
computational biology
```

grep’s default output often doesn’t give us enough context of a match when we need to inspect results by eye; only the matching line is printed to standard output. There are three useful options to get around this context before (-B), context: after (-A), and context before and after (-C). Each of these arguments takes how many lines of context to provide:

```
$ grep -B1 "AGATCGG" contam.fastq | head -n 6
@DJB775P1:248:D0MDGACXX:7:1202:12362:49613
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA
--
@DJB775P1:248:D0MDGACXX:7:1202:12782:49716
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
--
$ grep -A2 "AGATCGG" contam.fastq | head -n 6
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA
+
JJJJJIIJJJJJJHIHHHGHFFFFFFCEEEEEDBD?DDDDDDBDDDABDDCA
--
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
+
```

For example, if we wanted to find the Ensembl gene identifiers for both “Olfr1413” and “Olfr1411,” we
could use:

```
$ grep "Olfr141[13]" Mus_musculus.GRCm38.75_chr1_genes.txt
ENSMUSG00000058904 Olfr1413
ENSMUSG00000062497 Olfr1411
```

grep has an option to count how many lines match a pattern: -c. For example, suppose we wanted a quick look at how many genes start with “Olfr”:

```
$ grep -c "\tOlfr" Mus_musculus.GRCm38.75_chr1_genes.txt
27
```

For example, suppose we wanted to know how many small nuclear RNAs are in our Mus_musculus.GRCm38.75_chr1.gtf file.

```
$ grep -c 'gene_biotype "snRNA"' Mus_musculus.GRCm38.75_chr1.gtf
315
```

Suppose we wanted to extract all values of the “gene_id” field from the last column of our Mus_musculus.GRCm38.75_chr1.gtf file. This is easy with -o:

```
$ grep -E -o 'gene_id "\w+"' Mus_musculus.GRCm38.75_chr1.gtf | head -n 5
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000090025"
gene_id "ENSMUSG00000064842"
gene_id "ENSMUSG00000064842"
```

#### Cleaning a set of gene names with Unix data tools
```
grep -E -o 'gene_id "(\w+)"' Mus_musculus.GRCm38.75_chr1.gtf | \
cut -f2 -d" " | \
sed 's/"//g' | \
sort | \
uniq > mm_gene_id.txt
```

### Sorting Plain-Text Data with Sort

Very often we need to work with sorted plain-text data in bioinformatics. First, like cut, sort is designed to work with plain-text data with columns. Running sort without any arguments simply sorts a file alphanumerically by line:

```
$ cat example.bed
$ sort example.bed
```

However, using sort’s defaults of sorting alphanumerically by line doesn’t handle tabular data properly. So we need to sort by particular columns and to tell sort that certain columns are numeric values.

sort has a simple syntax to do this. Let’s look at how we’d sort example.bed by chromosome (first column), and start position (second column):

```
sort -k1,1 -k2,2n example.bed
chr1 9 28
chr1 10 19
chr1 26 39
chr1 32 47
chr1 40 49
chr2 35 54
chr3 11 28
chr3 16 27
```

We could then redirect the standard output stream of sort to a file:

```
$ sort -k1,1 -k2,2n example.bed > example_sorted.bed
```

Another example: 

```
$ sort -k1,1 -k4,4n Mus_musculus.GRCm38.75_chr1_random.gtf > \
Mus_musculus.GRCm38.75_chr1_sorted.gtf
```

We can check if a file is sorted according to our -k arguments using -c:

```
$ sort -k1,1 -k2,2n -c example_sorted.bed
$ echo $?
0
$ sort -k1,1 -k2,2n -c example.bed
sort: example.bed:4: disorder: chr1 40 49
$ echo $?
1
```

It’s also possible to sort in reverse order with the -r argument:

```
$ sort -k1,1 -k2,2n -r example.bed
chr3 11 28
chr3 16 27
chr2 35 54
chr1 9 28
chr1 10 19
chr1 26 39
chr1 32 47
chr1 40 49
```

If you’d like to only reverse the sorting order of a single column, you can append r on that column’s -k argument:

```
$ sort -k1,1 -k2,2nr example.bed
chr1 40 49
chr1 32 47
chr1 26 39
chr1 10 19
chr1 9 28
chr2 35 54
chr3 16 27
chr3 11 28
```

There are a few other useful sorting options to discuss, but these are available for GNU sort only (not the BSD version as found on OS X). The first is -V, which is a clever alphanumeric sorting routine that understands numbers inside strings.

```
$ sort -k1,1V -k2,2n example2.bed
chr1 34 49
chr10 30 42
chr10 31 47
chr11 6 16
chr2 15 19
chr2 17 22
chr2 27 46
chr22 32 46
```

Another option (only available in GNU sort) is to run sort with the --parallel option. For example, to use four cores to sort Mus_musculus.GRCm38.75_chr1_random.gtf:

```
$ sort -k1,1 -k4,4n --parallel 4 Mus_musculus.GRCm38.75_chr1_random.gtf
```

### Finding Unique Values in Uniq

Unix’s uniq takes lines from a file or standard input stream, and outputs all lines with consecutive duplicates removed.

Let’s first see an example of its behavior:

```
$ cat letters.txt
A
A
B
C
B
C
C
C
$ uniq letters.txt
A
B
C
B
C
```

As you can see, uniq does not return the unique values letters.txt—it only removes consecutive duplicate lines (keeping one). If instead we did want to find all unique lines in a file, we would first sort all lines using sort so that all identical lines are grouped next to each other, and then run uniq. For example:

```
$ sort letters.txt | uniq
A
B
C
```

uniq also has a tremendously useful option that’s used very often in command-line data processing: -c. This option shows the counts of occurrences next to the unique lines. For example:

```
$ uniq -c letters.txt
2 A
1 B
1 C
1 B
3 C
$ sort letters.txt | uniq -c
2 A
2 B
4 C
```

Both sort | uniq and sort | uniq -c are frequently used shell idioms in bioinformatics and worth memorizing. Combined with other Unix tools like grep and cut, sort and uniq can be used to summarize columns of tabular data:

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f3 | sort | uniq -c
25901 CDS
7588 UTR
36128 exon
2027 gene
2290 start_codon
2299 stop_codon
4993 transcript
```

If we wanted these counts in order from most frequent to least, we could pipe these results to sort -rn:

```
$ grep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf | cut -f3 | sort | uniq -c | \
sort -rn
36128 exon
25901 CDS
7588 UTR
4993 transcript
2299 stop_codon
2290 start_codon
2027 gene
```

If you want to see the number of features belonging to a particular gene identifier:

```
$ grep "ENSMUSG00000033793" Mus_musculus.GRCm38.75_chr1.gtf | cut -f3 | sort \
| uniq -c
13 CDS
3 UTR
14 exon
1 gene
1 start_codon
1 stop_codon
1 transcript
```

uniq can also be used to check for duplicates with the -d option. With the -d option, uniq outputs duplicated lines only.

```
$ uniq -d mm_gene_names.txt
# no output
$ uniq -d mm_gene_names.txt | wc -l
0
```

A file with duplicates, like the test.bed file, has multiple lines returned:

```
$ uniq -d test.bed | wc -l
22925
```

### Join

The Unix tool join is used to join different files together by a common column. This is easiest to understand with simple test data. Let’s use our example.bed BED file, and example_lengths.txt, a file containing the same chromosomes as example.bed with their lengths. Both files look like this:

```
$ cat example.bed
chr1 26 39
chr1 32 47
chr3 11 28
chr1 40 49
chr3 16 27
chr1 9 28
chr2 35 54
chr1 10 19
$ cat example_lengths.txt
chr1 58352
chr2 39521
chr3 24859
```

Our goal is to append the chromosome length alongside each feature (note that the result will not be a valid BED-formatted file, just a tab-delimited file).

To append the chromosome lengths to example.bed, we first need to sort both files by the column to be joined on.

```
$ sort -k1,1 example.bed > example_sorted.bed
$ sort -c -k1,1 example_lengths.txt # verifies is already sorted
```

Then, to sort them, we can use:

```
$ join -1 1 -2 1 example_sorted.bed example_lengths.txt > example_with_lengths.txt
$ cat example_with_lengths.txt
chr1 10 19 58352
chr1 26 39 58352
chr1 32 47 58352
chr1 40 49 58352
chr1 9 28 58352
chr2 35 54 39521
chr3 11 28 24859
chr3 16 27 24859
```

But look what happens if our second file, example_lengths.txt, is truncated such
that it doesn’t have the lengths for chr3:

```$ head -n2 example_lengths.txt > example_lengths_alt.txt # truncate file
$ join -1 1 -2 1 example_sorted.bed example_lengths_alt.txt
chr1 10 19 58352
chr1 26 39 58352
chr1 32 47 58352
chr1 40 49 58352
chr1 9 28 58352
chr2 35 54 39521
$ join -1 1 -2 1 example_sorted.bed example_lengths_alt.txt | wc -l
6
```

GNU join implements the -a option to include unpairable lines—ones that do not have an entry in either file.

```
$ join -1 1 -2 1 -a 1 example_sorted.bed example_lengths_alt.txt # GNU join only
chr1 10 19 58352
chr1 26 39 58352
chr1 32 47 58352
chr1 40 49 58352
chr1 9 28 58352
chr2 35 54 39521
chr3 11 28
chr3 16 27
```

### Text Processing with Awk

We’ll introduce the basics of Awk in this section—enough to get you started with using Awk in bioinformatics.

First, Awk processes input data a record at a time. Each record is composed of fields, separate chunks that Awk automatically separates. Because Awk was designed to work with tabular data, each record is a line, and each field is a column’s entry for that record.

First, we can simply mimic cat by omitting a pattern and printing an entire record with the variable $0:

```
$ awk '{ print $0 }' example.bed
chr1 26 39
chr1 32 47
chr3 11 28
chr1 40 49
chr3 16 27
chr1 9 28
chr2 35 54
chr1 10 19
```

Awk supports arithmetic with the standard operators +, -, *, /, % (remainder), and ^ (exponentiation). We can subtract within a pattern to calculate the length of a feature, and filter on that expression:

```
$ awk '$3 - $2 > 18' example.bed
chr1 9 28
chr2 35 54
```

We can also chain patterns, by using logical operators && (AND), || (OR), and ! (NOT). For example, if we wanted all lines on chromosome 1 with a length greater than 10:

```
$ awk '$1 ~ /chr1/ && $3 - $2 > 10' example.bed
chr1 26 39
chr1 32 47
chr1 9 28
```

We would have to take the sum feature lengths, and then divide by the total number of records. We can do this with:

```
$ awk 'BEGIN{ s = 0 }; { s += ($3-$2) }; END{ print "mean: " s/NR };' example.bed
mean: 14
```

We can use NR to extract ranges of lines, too; for example, if we wanted to extract all lines between 3 and 5 (inclusive):

```
$ awk 'NR >= 3 && NR <= 5' example.bed
chr3 11 28
chr1 40 49
chr3 16 27
```
We could generate a three-column BED file from Mus_musculus.GRCm38.75_chr1.gtf as follows:

```
$ awk '!/^#/ { print $1 "\t" $4-1 "\t" $5 }' Mus_musculus.GRCm38.75_chr1.gtf | \
head -n 3
1 3054232 3054733
1 3054232 3054733
1 3054232 3054733
```

Awk also has a very useful data structure known as an associative array. Associative arrays behave like Python’s dictionaries or hashes in other languages. We can create an associative array by simply assigning a value to a key. For example, suppose we
wanted to count the number of features (third column belonging to the gene “Lypla1.” We could do this by incrementing their values in an associative array:

```
$ awk '/Lypla1/ { feature[$3] += 1 }; \
END { for (k in feature) \
print k "\t" feature[k] }' Mus_musculus.GRCm38.75_chr1.gtf
```

It’s worth noting that there’s an entirely Unix way to count features of a particular gene: grep, cut, sort, and uniq -c:

```
$ grep "Lypla1" Mus_musculus.GRCm38.75_chr1.gtf | cut -f 3 | sort | uniq -c
56 CDS
24 UTR
69 exon
1 gene
5 start_codon
5 stop_codon
9 transcript
```

### Bioawk: An Awk for Biological Formats

Imagine extending Awk’s powerful processing of tabular data to processing tasks involving common bioinformatics formats like FASTA/FASTQ, GTF/GFF, BED, SAM, and VCF. Let’s look at Bioawk’s supported
input formats and what variables these formats set:

To install bioawk on Ubuntu we can follow this commands: https://silico-sciences.com/2015/12/install-bioawk-on-ubuntu/

```
$ bioawk -c help
```

Bioawk is also quite useful for processing FASTA/FASTQ files. For example, we could use it to turn a FASTQ file into a FASTA file:

```
$ bioawk -c fastx '{print ">"$name"\n"$seq}' contam.fastq | head -n 4
>DJB775P1:248:D0MDGACXX:7:1202:12362:49613
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA
>DJB775P1:248:D0MDGACXX:7:1202:12782:49716
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
```

Bioawk can also serve as a method of counting the number of FASTQ/FASTA entries:

```
$ bioawk -c fastx 'END{print NR}' contam.fastq
8
```

Bioawk is also useful for creating a table of sequence lengths from a FASTA file. For example, to create a table of all chromosome lengths of the Mus musculus genome:

```
$ bioawk -c fastx '{print $name,length($seq)}' \
Mus_musculus.GRCm38.75.dna_rm.toplevel.fa.gz > mm_genome.txt
$ head -n 4 mm_genome.txt
1 195471971
10 130694993
11 122082543
12 120129022
13 120421639
14 124902244
```

If we wanted to return all variants for which individuals ind_A and ind_B have identical genotypes (note that this assumes a fixed allele order like ref/alt or major/minor):

```
$ bioawk -c hdr '$ind_A == $ind_B {print $id}' genotypes.txt
S_001
S_003
S_005
S_008
S_009
```

### Advanced Shell Tricks

Now we’re ready to dig into a few more advanced shell tricks.

#### Subshells

The first trick we’ll cover is using Unix subshells. Before explaining this trick, it’s helpful to remember the difference between sequential commands (connected with && or ;), and piped commands (connected with |). Sequential commands are simply run one after the other; the previous command’s standard output does not get passed to the next program’s standard in. In contrast, connecting two programs with pipes means the first program’s standard out will be piped into the next program’s
standard in. If we use command1 && command2, command2 will only run if command1 completed with a zero-exit status.

Subshells allow us to execute sequential com‐
mands together in a separate shell process. This is useful primarily to group sequential commands together (such that their output is a single stream. This gives us a new way to construct clever one-liners and has practical uses in command-line data processing. Let’s look at a toy example first:

```
$ echo "this command"; echo "that command" | sed 's/command/step/'
this command
that step
$ (echo "this command"; echo "that command") | sed 's/command/step/'
this step
that step
```

Consider the problem of sorting a GTF file with a metadata header. We can’t simply sort the entire file with sort, because this header could get shuffled in with rows of data. Instead, we want to sort everything except the header, but still include the header at the top of the final sorted file. We can solve this problem using a subshell to group sequential commands that print the header to standard out and sort all other lines by chromosome and start position, printing all lines to standard out after the header.

```
$ (zgrep "^#" Mus_musculus.GRCm38.75_chr1.gtf.gz; \
zgrep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf.gz | \
sort -k1,1 -k4,4n) | less
```

Because we’ve used a subshell, all standard output from these sequential commands will be combined into a single stream, which here is piped to less. To write this stream to a file, we could redirect this stream to a file using something like > 
Mus_musculus.GRCm38.75_chr1_sorted.gtf. But a better approach would be to use gzip to compress this stream before writing it to disk:

```
$ (zgrep "^#" Mus_musculus.GRCm38.75_chr1.gtf.gz; \
zgrep -v "^#" Mus_musculus.GRCm38.75_chr1.gtf.gz | \
sort -k1,1 -k4,4n) | gzip > Mus_musculus.GRCm38.75_chr1_sorted.gtf.gz
```

### Named Pipes and Process Substitution

Throughout this chapter, we’ve used pipes to connect command-line tools to build custom data-processing pipelines. However, some programs won’t interface with the Unix pipes we’ve come to love and depend on. For example, certain bioinformatics tools read in multiple input files and write to multiple output files:

```
$ processing_tool --in1 in1.fq --in2 in2.fq --out1 out2.fq --out2.fq
```

In this case, the imaginary program processing_tool requires two separate input files, and produces two separate output files. Because each file needs to be provided separately, we can’t pipe the previous processing step’s results through processing_tool’s standard in.

A named pipe, also known as a FIFO (First In First Out, a concept in computer science), is a special sort of file. Regular pipes are anonymous—they don’t have a name, and only persist while both processes are running. Named pipes behave like files, and are persistent on your filesystem. We can create a named pipe with the program mkfifo:

```
$ mkfifo fqin
$ ls -l fqin
prw-r--r-- 1 vinceb staff 0 Aug 5 22:50 fqin
```

As a toy example, we can simulate this by using echo to redirect some text into a named pipe (running it in the background, so we can have our prompt back), and then cat to read the data back out:

```
$ echo "hello, named pipes" > fqin &
[1] 16430
$ cat fqin
[1] + 16430 done
hello, named pipes
$ rm fqin
```

Although the syntax is similar to shell redirection to a file, we’re not actually writing anything to our disk. Named pipes provide all of the computational benefits of pipes with the flexibility of interfacing with files.

However, creating and removing these file-like named pipes is a bit tedious. Programmers like syntactic shortcuts, so there’s a way to use named pipes without having to explicitly create them. This is called process substitution, or sometimes known as
anonymous named pipes. These allow you to invoke a process, and have its standard output go directly to a named pipe.

If we were to re-create the previous toy example with process substitution, it would look as follows:

```
$ cat <(echo "hello, process substitution")
hello, process substitution
```

Order to follow: chapter 1, 2, 3, 6, 7, 9, 10, 11, 12, 8.

## Chapter 8 - A Rapid Introduction to the R Language

(175 to 261)

To install R in all OSs, we can follow this tutorial: 
https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu#ubuntu

## Chapter 9 - Working with Range Data

Luckily for bioinformaticians, every genome from every branch of life on earth consists of chromosome sequences that can be represented on a computer in the same way: as a set of nucleotide sequences (genomic variation and assembly uncertainty
aside). Each separate sequence represents a reference DNA molecule, which may correspond to a fully assembled chromosome, or a scaffold or contig in a partially assembled genome. Although nucleotide sequences are linear, they may also represent biologically circular chromosomes (e.g., with plasmids or mitochondria) that have been cut.

Many types of genomic data are linked to a specific genomic region, and this region can be represented as a range containing consecutive positions on a chromosome. Annotation data and genomic features like gene models, genetic variants like SNPs and insertions/deletions, transposable elements, binding sites, and statistics like pairwise diversity and GC content can all be represented as ranges on a linear chromosome sequence. Sequencing read alignment data resulting from experiments like whole genome resequencing, RNA-Seq, ChIP-Seq, and bisulfite sequencing can also be represented as ranges.

Once our genomic data is represented as ranges on chromosomes, there are numerous range operations at our disposal to tackle tasks like finding and counting overlaps, calculating coverage, finding nearest ranges, and extracting nucleotide sequences from specific ranges.

### A Crash Course in Genomic Ranges and Coordinate Systems

So what are ranges exactly? Ranges are integer intervals that represent a subsequence of consecutive positions on a sequence like a chromosome.

Ranges alone only specify a region along a single sequence like a chromosome; to specify a genomic region or position, we need three necessary pieces of
information:

- **Chromosome name**: This is also known as sequence name (to allow for sequences that aren’t fully assembled, such as scaffolds or contigs). Each

- **Range**: For example, 114,414,997 to 114,693,772 or 3,173,498 to 3,179,449. Ranges are how we specify a single subsequence on a chromosome sequence.

- **Strand**: Because chromosomal DNA is double-stranded, features can reside on either the forward (positive) or reverse (negative) strand. Many features on a chromosome are strand-specific, so we need to specify which strand these features are on.

Despite the convenience that comes with representing and working with genomic ranges, there are unfortunately some gritty details we need to be aware of. First, there are two different flavors of range systems used by bioinformatics data formats and software programs:

- **0-based coordinate system, with half-closed, half-open intervals**: the first base of a sequence is position 0 and the last base’s position is the length of the sequence - 1 (exclude the last position). Example: [1, 5). In fact, Pytho strings work like this to:

```
$ python3
>>> "CTTACTTCGAAGGCTG"[1:5]
'TTAC'
```

- **1-based coordinate system, with closed intervals**: As you might have guessed, with 1-based systems the first base of a sequence is given the position 1. Because positions are counted as we do natural numbers, the last position in a sequence is always equal to its length. Closed intervals means like: [2,5]. 

R uses 1-based indexing for its vectors and strings, and extracting a portion of a string with substr() uses closed intervals:

```
> substr("CTTACTTCGAAGGCTG", 2, 5)
[1] "TTAC"
```

Here, we can see the range types of common bioinformatics formats:

| Format/library  | Type |
| ------------- | ------------- |
| BED  | 0-based  |
| GTF  | 1-based  |
| GFF  | 1-based  |
| SAM  | 1-based  |
| BAM  | 0-based  |
| VCF  | 1-based  |
| BCF  | 0-based  |
| Wiggle  | 1-based  |
| GenomicRanges  | 1-based  |
| BLAST  | 1-based  |
| GenBank/EMBL Feature Table  | 1-based  |

The second gritty detail we need to worry about is strand. There’s little to say except: you need to mind strand in your work. Because DNA is double stranded, genomic features can lie on either strand.
However, a genomic feature can be either on the forward or reverse strand. For genomic features like protein coding regions, strand matters and must be specified. For example, a range representing a protein coding region only makes biological sense given the appropriate strand.

### An Interactive Introduction to Range Data with GenomicRanges

(Continue later)

## Chapter 10 - Working with Sequence Data

Nucleotide (and protein) sequences are stored in two plain-text formats widespread in bioinformatics: FASTA and FASTQ. We’ll discuss each format and their limitations in this section, and then see some tools for working with data in these formats.

### The FASTA Format

The FASTA format originates from the FASTA alignment suite, created by William R. Pearson and David J. Lipman. The FASTA format is used to store any sort of sequence data not requiring per-base pair quality scores. This includes reference genome files, protein sequences, coding DNA sequences (CDS), transcript sequences, and so on.

FASTA files are composed of sequence entries, each containing two parts: a description and the sequence data. The description line begins with a greater than symbol (>) and contains the sequence identifier and other (optional) information. The sequence data begins on the next line after the description, and continues until there’s another description line (a line beginning with >) or the file ends.

The egfr_flank.fasta file in this chapter’s GitHub directory is an example FASTA file:

```
$ head -10 egfr_flank.fasta
>ENSMUSG00000020122|ENSMUST00000138518
CCCTCCTATCATGCTGTCAGTGTATCTCTAAATAGCACTCTCAACCCCCGTGAACTTGGT
TATTAAAAACATGCCCAAAGTCTGGGAGCCAGGGCTGCAGGGAAATACCACAGCCTCAGT
TCATCAAAACAGTTCATTGCCCAAAATGTTCTCAGCTGCAGCTTTCATGAGGTAACTCCA
GGGCCCACCTGTTCTCTGGT
>ENSMUSG00000020122|ENSMUST00000125984
GAGTCAGGTTGAAGCTGCCCTGAACACTACAGAGAAGAGAGGCCTTGGTGTCCTGTTGTC
TCCAGAACCCCAATATGTCTTGTGAAGGGCACACAACCCCTCAAAGGGGTGTCACTTCTT
CTGATCACTTTTGTTACTGTTTACTAACTGATCCTATGAATCACTGTGTCTTCTCAGAGG
CCGTGAACCACGTCTGCAAT
```

The FASTA format’s simplicity and flexibility comes with an unfortunate downside: the FASTA format is a loosely defined ad hoc format (which unfortunately are quite common in bioinformatics). This is why it’s usually preferable to use existing FASTA/FASTQ parsing libraries instead of implementing your own; existing libraries have already been vetted by the open source community.

Most troubling about the FASTA format is that there’s no universal specification for the format of an identifier in the description. For example, should the following FASTA descriptions refer to the same entry?

```
>ENSMUSG00000020122|ENSMUST00000138518
> ENSMUSG00000020122|ENSMUST00000125984
>ENSMUSG00000020122|ENSMUST00000125984|epidermal growth factor receptor
>ENSMUSG00000020122|ENSMUST00000125984|Egfr
>ENSMUSG00000020122|ENSMUST00000125984|11|ENSFM00410000138465
```

Without a standard scheme for identifiers, we can’t use simple exact matching to check if an identifier matches a FASTA entry header line. Instead, we would need to rely on fuzzy matching between FASTA descriptions and our identifier.

A common naming convention is to split the description line into two parts at the first space: the identifier and the comment. A sequence in this format would look like:

```
>gene_00284728 length=231;type=dna
GAGAACTGATTCTGTTACCGCAGGGCATTCGGATGTGCTAAGGTAGTAATCCATTATAAGTAACATGCGCGGAATATCCG
GAGGTCATAGTCGTAATGCATAATTATTCCCTCCCTCAGAAGGACTCCCTTGCGAGACGCCAATACCAAAGACTTTCGTA
GCTGGAACGATTGGACGGCCCAACCGGGGGGAGTCGGCTATACGTCTGATTGCTACGCCTGGACTTCTCTT
```

### The FASTQ Format

The FASTQ format extends FASTA by including a numeric quality score to each base in the sequence. The FASTQ format is widely used to store high-throughput sequencing data, which is reported with a per-base quality score indicating the confidence of each base call.

The FASTQ format looks like:

```
@DJB775P1:248:D0MDGACXX:7:1202:12362:49613 (1)
TGCTTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAA (2)
+ (3)
JJJJJIIJJJJJJHIHHHGHFFFFFFCEEEEEDBD?DDDDDDBDDDABDDCA (4)
@DJB775P1:248:D0MDGACXX:7:1202:12782:49716
CTCTGCGTTGATACCACTGCTTACTCTGCGTTGATACCACTGCTTAGATCGG
+
IIIIIIIIIIIIIIIHHHHHHFFFFFFEECCCCBCECCCCCCCCCCCCCCCC
```

1. Description line, which contains the identifier and other information
2. Sequence data.
3. End of the sequence.
4. Quality data encoded with ASCII characters.

We can use the fact that the number of quality score characters must be equal to the number of sequence characters to reliably parse this format—which is how the readfq parser introduced later on works, because @ character is not only the delimiter of the header but can be a quality character.

#### Counting FASTA/FASTQ Entries

```
$ grep -c "^>" egfr_flank.fasta
5
```

This approach will fail on FASTQ, because it counts que @ quality character too. 

```
$ grep -c "^@" untreated1_chr4.fq
208779
```

If you’re absolutely positive your FASTQ file uses four lines per sequence entry, you can estimate the number of sequences by estimating the number of lines with wc -l and dividing by four. If you’re unsure if some of your FASTQ entries wrap across many lines, a more robust way to count sequences is with bioawk:

```
$ bioawk -cfastx 'END{print NR}' untreated1_chr4.fq
204355
```

### Nucleotide Codes

With the basic FASTA/FASTQ formats covered, let’s look at the standards for encoding nucleotides and base quality scores in these formats. Clearly, encoding nucleotides is simple: A, T, C, G. Lowercase bases are often used to indicate soft masked repeats or low complexity sequences.
Repeats and low-complexity sequences may also be hard masked, where nucleotides are replaced with N (or sometimes an X). Degenerate (or ambiguous) nucleotide codes are used to represent two or more bases. For example, N is used to represent any base.

### Base Qualities

Each sequence base of a FASTQ entry has a corresponding numeric quality score in the quality line(s). Each base quality scores is encoded as a single ASCII character. Quality lines look like a string of random characters, like the fourth line here:

```
@AZ1:233:B390NACCC:2:1203:7689:2153
GTTGTTCTTGATGAGCCATGAGGAAGGCATGCCAAATTAAAATACTGGTGCGAATTTAAT
+
CCFFFFHHHHHJJJJJEIFJIJIJJJIJIJJJJCDGHIIIGIGIJIJIIIIJIJJIJIIH
```

Remember, ASCII characters are just represented as integers between 0 and 127 under the hood. Because not all ASCII characters are printable to screen (e.g., character echoing "\07" makes a “ding” noise)qualities are restricted to the printable ASCII characters, ranging from 33 to 126 (the space character, 32, is omitted).

All programming languages have functions to convert from a character to its decimal ASCII representation, and from ASCII decimal to character. In Python, these are the functions ord() and chr(), respectively. Let’s use ord() in Python’s interactive interpreter to convert the quality characters to a list of ASCII decimal representations:

```
>>> qual = "JJJJJJJJJJJJGJJJJJIIJJJJJIGJJJJJIJJJJJJJIJIJJJJHHHHHFFFDFCCC"
>>> [ord(b) for b in qual]
[74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 71, 74, 74, 74, 74, 74, 73,
73, 74, 74, 74, 74, 74, 73, 71, 74, 74, 74, 74, 74, 73, 74, 74, 74, 74, 74,
74, 74, 73, 74, 73, 74, 74, 74, 74, 72, 72, 72, 72, 72, 70, 70, 70, 68, 70,
67, 67, 67]
```

Unfortunately, converting these ASCII values to meaningful quality scores can be tricky because there are three different quality schemes: Sanger, Solexa, and Illumina. The Open Bioinformatics Foundation (OBF), which is responsible for projects like Biopython, BioPerl, and BioRuby, gives these the names fastq-sanger, fastq-solexa, and fastq-illumina.

Below, we can see FASTQ quality schemes:

| Name  | ASCII range | Offset | Quality score type | Quality score range
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Sanger, Illumina (versions 1.8 onward) | 33–126 | 33 | PHRED | 0–93
| Solexa, early Illumina (before 1.3) | 59–126  | 64 | Solexa | 5-62 
| Illumina (versions 1.3–1.7) | 64–126  | 64 | PHRED | 0-62

First, we need to subtract an offset to convert this Sanger quality score to a PHRED quality score. PHRED was an early base caller written by Phil Green, used for fluorescence trace data written by Phil Green. Looking at the table above, notice that the Sanger
format’s offset is 33, so we subtract 33 from each quality score:

```
>>> phred = [ord(b)-33 for b in qual]
>>> phred
[41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 38, 41, 41, 41, 41, 41, 40,
40, 41, 41, 41, 41, 41, 40, 38, 41, 41, 41, 41, 41, 40, 41, 41, 41, 41, 41,
41, 41, 40, 41, 40, 41, 41, 41, 41, 39, 39, 39, 39, 39, 37, 37, 37, 35, 37,
34, 34, 34]
```

Using the formula to convert quality score to the estimated probability the base is correct:

```
>>> [10**(-q/10) for q in phred]
[1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05,
1e-05, 0.0001, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 0.0001, 0.0001, 1e-05,
1e-05, 1e-05, 1e-05, 1e-05, 0.0001, 0.0001, 1e-05, 1e-05, 1e-05, 1e-05,
1e-05, 0.0001, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 0.0001,
1e-05, 0.0001, 1e-05, 1e-05, 1e-05, 1e-05, 0.0001, 0.0001, 0.0001, 0.0001,
0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]
```

### Example: Inspecting and Trimming Low-Quality Bases

Our Python list of base accuracies is useful as a learning tool to see how to convert qualities to probabilities, but it won’t help us much to understand the quality profiles of millions of sequences. In this sense, a picture is worth a thousand words—and there’s software to help us see the quality distributions across bases in reads. We’ll use qrqc R package through Bioconductor in examples so we can tinker with how we visualize this data ourselves.

Let’s first install all the necessary programs for this example. First, install qrqc in R with:

```
> library(BiocInstaller)
> biocLite('qrqc')
```

Next, let’s install two programs that will allow us to trim low-quality bases: **sickle** and **seqtk**. seqtk is a general-purpose sequence toolkit written by Heng Li that contains a subcommand for trimming low-quality bases off the end of sequences.

```
sudo apt-get install sickle
sudo apt-get install seqtk
```

After getting these programs installed, let’s trim the untreated1_chr4.fq FASTQ file in this chapter’s directory in the GitHub repository. This FASTQ file was generated from the untreated1_chr4.bam BAM file in the pasillaBamSubset Bioconductor package. To keep
things simple, we’ll use each program’s default settings. Starting with sickle:

```
$ sickle se -f untreated1_chr4.fq -t sanger -o untreated1_chr4_sickle.fq
FastQ records kept: 203121
FastQ records discarded: 1234
```

Now, let’s run seqtk trimfq, which takes a single argument and outputs trimmed sequences through standard out:

$ seqtk trimfq untreated1_chr4.fq > untreated1_chr4_trimfq.fq

Let’s compare these results in R. We’ll use qrqc to collect the distributions of quality by position in these files, and then visualize these using ggplot2. We could load these in one at a time, but a nice workflow is to automate this with lapply(). The script is on chapter10 package and produces two plots (plot1 and plot2): where We see the effect both trimming programs have on our data’s quality distributions in  low-quality bases, we narrow the quality ranges in base positions further in the read.
Furthermore, we see this increases mean quality across across the read, but we still
see a decline in base quality over the length of the reads.

In one line, we can trim low-quality bases from the ends of these sequences—running the trimming commands is not difficult. The more important step is to visualize what these trimming programs did to our data by comparing the files before and after
trimming.

### A FASTA/FASTQ Parsing Example: Counting Nucleotides

It’s not too difficult to write your own FASTA/FASTQ parser and it’s a useful educational programming exercise. But when it comes to using a parser for processing real data, it’s best to use a reliable existing library.

We’ll use Heng Li’s readfq implementation because it parses both FASTA and FASTQ files, it’s simple to use, and is standalone (meaning it doesn’t require installing any dependencies). **Biopython** and **BioPerl** are two popular libraries with good alternative FASTA/FASTQ parsers.

We’ll use the Python implementation of readfq, readfq.py to read FASTQ files and let’s write a simple program that counts the number of each IUPAC nucleotide in a file (nuccont.py):

This version takes input through standard in, so after saving this file and adding execute permissions with **chmod +x nuccount.py**, we could run it with:

```
$ cat contam.fastq | ./nuccount.py
A 103
C 110
G 94
T 109
R 0
Y 0
S 0
W 0
K 0
M 0
B 0
D 0
H 0
V 0
N 0
```

There are many improvements we could add to this script: add support for persequence base composition statistics, take file arguments rather than input through standard in, count soft-masked (lowercase) characters, count CpG sites, or warn when a non-IUPAC nucleotide is found in a file. Counting nucleotides is simple—the most complex part of the script is readfq(). This is the beauty of reusing code: well-written functions and libraries prevent us from having to rewrite complex parsers.

## Chapter 11 - Working with Alignment Data

In Chapter 9, we learned about range formats such as BED and GTF, which are often used to store genomic range data associated with genomic feature annotation data such as gene models. Other kinds of range-based formats are designed for storing large amounts of alignment data—for example, the results of aligning millions (or billions) of high-throughput sequencing reads to a genome. In this chapter, we’ll look at the most common high-throughput data alignment format: the Sequence Alignment/Mapping (SAM) format for mapping data (and its binary analog, BAM). The SAM and BAM formats are the standard formats for storing sequencing reads mapped to a reference.

We study SAM and BAM for two reasons. First, a huge part of bioinformatics work is manipulating alignment files. Nearly every high-throughput sequencing experiment involves an alignment step that produces alignment data in the SAM/BAM formats.
Second, the skills developed through learning to work with SAM/BAM files are extensible and more widely applicable than to these specific formats.

### Getting to Know Alignment Formats: SAM and BAM

Before learning to work with SAM/BAM, we need to understand the structure of these formats. We’ll do this by using celegans.sam, a small example SAM file.

#### SAM header

Files in the SAM format consist of a header section and an alignment section. Because SAM files are plain text (unlike their binary counterpart, BAM), we can take a peek at a few lines of the header with head:

```
$ head -n 10 celegans.sam
@SQ SN:I LN:15072434
@SQ SN:II LN:15279421
@SQ SN:III LN:13783801
@SQ SN:IV LN:17493829
@SQ SN:MtDNA LN:13794
@SQ SN:V LN:20924180
@SQ SN:X LN:17718942
@RG ID:VB00023_L001 SM:celegans-01
@PG ID:bwa PN:bwa VN:0.7.10-r789 [...]
I_2011868_2012306_0:0:0_0:0:0_2489 83 I 2012257 40 50M [...]
```

The standard way of interacting with SAM/BAM files is through the SAMtools command-line program (samtools), which we’ll use extensively throughout the rest of this chapter. 

A universal way to look at an entire SAM/BAM header is with samtools view option -H:

```
$ sudo apt-get install samtools
$ samtools view -H celegans.sam
@SQ SN:I LN:15072434
@SQ SN:II LN:15279421
@SQ SN:III LN:13783801
[...]
```

This also works with BAM files:

```
$ samtools view -H celegans.bam
@SQ SN:I LN:15072434
@SQ SN:II LN:15279421
@SQ SN:III LN:13783801
[...]
```

Of course, all our usual Unix tricks can be combined with samtools commands through piping results to other commands. For example, we could see all read groups with:

```
$ samtools view -H celegans.bam | grep "^@RG"
@RG ID:VB00023_L001 SM:celegans-01
```

#### SAM alignment section

The alignment section contains read alignments (and usually includes reads that did not align, but this depends on the aligner and file). Each alignment entry is composed of 11 required fields (and optional fields after this).

```
$ samtools view celegans.sam | tr '\t' '\n' | head -n 11
I_2011868_2012306_0:0:0_0:0:0_2489
83
I
2012257
40
50M
=
2011868
-439
CAAAAAATTTTGAAAAAAAAAATTGAATAAAAATTCACGGATTTCTGGCT
22222222222222222222222222222222222222222222222222
```

#### CIGAR strings

Like bitwise flags, SAM’s CIGAR strings are another specialized way to encode information about an aligned sequence. While bitwise flags store true/false properties about an alignment, CIGAR strings encode information about which bases of an alignment are matches/mismatches, insertions, deletions, soft or hard clipped, and so on. I’ll assume you are familiar with the idea of matches, mismatches, insertions, and deletions, but it’s worth describing soft and hard clipping (as SAM uses them).

The cigar operations are:

M 0 Alignment match (note that this could
be a sequence match or mismatch!)
I 1 Insertion (to reference)
D 2 Deletion (from reference)
N 3 Skipped region (from reference)
S 4 Soft-clipped region (soft-clipped
regions are present in sequence in SEQ field)
H 5 Hard-clipped region (not in sequence in SEQ field)
P 6 Padding (see section 3.1 of the SAM
format specication for detail)
= 7 Sequence match

For example, a fully aligned 51 base pair read without insertions or deletions would
have a CIGAR string containing a single length/operation pair: 51M.

#### Mapping Qualities

Our discussion of the SAM and BAM formats is not complete without mentioning mapping qualities (Li et al., 2008). Mapping qualities are one of the most important diagnostics in alignment. All steps downstream of alignment in all bioinformatics projects (e.g., SNP calling and genotyping, RNA-seq, etc.) critically depend on reliable mapping. Mapping qualities quantify mapping reliability by estimating how likely a read is to actually originate from the position the aligner has mapped it to.

We can use mapping qualities to filter out likely incorrect alignments (which we can do with samtools view, which we’ll learn about later), find regions where mapping quality is unusually low among most alignments (perhaps in repetitive or paralogous
regions), or assess genome-wide mapping quality distributions (which could indicate alignment problems in highly repetitive or polyploid genomes).

### Command-Line Tools for Working with Alignments in the SAM Format

In this section, we’ll learn about the Samtools suite of tools for manipulating and working with SAM,BAM, and CRAM files. These tools are incredibly powerful, and becoming skilled in working with these tools will allow you to both quickly move forward in file-processing tasks and explore the data in alignment files.

#### Using samtools view to Convert between SAM and BAM

Many samtools subcommands such as sort, index, depth,and mpileup all require input files (or streams) to be in BAM format for efficiency, so we often need to convert between plain-text SAM and binary BAM formats. samtools view allows us to convert SAM to BAM with the -b option:

```
$ samtools view -b celegans.sam > celegans_copy.bam
```

Similarly, we can go from BAM to SAM:

```
$ samtools view celegans.bam > celegans_copy.sam
$ head -n 3 celegans_copy.sam
I_2011868_2012306_0:0:0_0:0:0_2489 83 I 2012257 40 [...]
I_2011868_2012306_0:0:0_0:0:0_2489 163 I 2011868 60 [...]
I_13330604_13331055_2:0:0_0:0:0_3dd5 83 I 13331006 60 [...]
```

However, samtools view will not include the SAM header (see “The SAM Header” on page 356) by default. SAM files without headers cannot be turned back into BAM files:

```
$ samtools view -b celegans_copy.sam > celegans_copy.bam
[E::sam_parse1] missing SAM header
[W::sam_read1] parse error at line 1
[main_samview] truncated file.
```

Converting BAM to SAM loses information when we don’t include the header. We can include the header with -h:

```
$ samtools view -h celegans.bam > celegans_copy.sam
$ samtools view -b celegans_copy.sam > celegans_copy.bam #now we can convert back
```

#### Samtools Sort and Index

We sort alignments by their alignment position with samtools sort:

```
$ samtools sort celegans_unsorted.bam celegans_sorted
```

We can use the samtools sort option -m to increase the memory, and -@ to specify how many threads to use. For example:

```
$ samtools sort -m 4G -@ 2 celegans_unsorted.bam celegans_sorted
```

Often, we want to work with alignments within a particular region in the genome. For example, we may want to extract these reads using samtools view or only call SNPs within this region. Iterating through an entire BAM file just to work with a subset of reads at a position would be inefficient; consequently, BAM files can be indexed. The BAM file must be sorted first, and we cannot index SAM files. To index a position-sorted BAM file, we simply use:

```
$ samtools index celegans_sorted.bam
```

(Continue later: 368 to 394)

## Chapter 12 - Bioinformatics Shell Scripting, Writing Pipelines, and Parallelizing Tasks

I’ve waited until the penultimate chapter this book to share a regrettable fact: everyday bioinformatics work often involves a great deal of tedious data processing. Bioinformaticians regularly need to run a sequence of commands on not just one file, but
dozens (sometimes even hundreds) of files. Consequently, a large part of bioinformatics is patching together various processing steps into a pipeline, and then repeatedly applying this pipeline to many files. This isn’t exciting scientific work, but it’s a necessary hurdle before tackling more exciting analyses.

While writing pipelines is a daily burden of bioinformaticians, it’s essential that pipelines are written to be robust and reproducible. Pipelines must be robust to problems that might occur during data processing. When we execute a series of commands on data directly into the shell, we usually clearly see if something goes awry—output files
are empty when they should contain data or programs exit with an error. However, when we run data through a processing pipeline, we sacrifice the careful attention we paid to each step’s output to gain the ability to automate processing of numerous files. The catch is that not only are errors likely to still occur, they’re more likely to occur because we’re automating processing over more data files and using more steps. For these reasons, it’s critical to construct robust pipelines.

Similarly, pipelines also play an important role in reproducibility. A well-crafted pipeline can be a perfect record of exactly how data was processed. In the best cases, an individual could download your processing scripts and data, and easily replicate your exact steps.

In this chapter, we’ll learn the essential tools and skills to construct robust and repro‐

ducible pipelines. We’ll see how to write rerunnable Bash shell scripts, automate file-processing tasks with find and xargs, run pipelines in parallel, and see a simple makefile.

### Basic Bash Scripting

Bash, the shell we’ve used interactively throughout the book, is also a full-fledged scripting language. Like many other tools presented in this book, the trick to using Bash scripts effectively in bioinformatics is knowing when to use them and when not to. Unlike Python, Bash is not a general-purpose language. Bash is explicitly designed to make running and interfacing command-line programs as simple as possible (a good characteristic of a shell!). For these reasons, Bash often takes the role as the duct tape language of bioinformatics (also referred to as a glue language), as it’s used to tape many commands together into a cohesive workflow.

Before digging into how to create pipelines in Bash, it’s important to note that Python may be a more suitable language for commonly reused or advanced pipelines. Python is a more modern, fully featured scripting language than Bash. Compared to Python,
Bash lacks several nice features useful for data-processing scripts: better numeric type
support, useful data structures, better string processing, refined option parsing, availability of a large number of libraries, and powerful functions that help with structuring your programs. However, there’s more overhead when calling command-line
programs from a Python script (known as calling out or shelling out) compared to Bash. Although Bash lacks some of Python’s features, Bash is often the best and quickest “duct tape” solution (which we often need in bioinformatics).

#### Writing and Running Robust Bash Scripts

Most Bash scripts in bioinformatics are simply commands organized into a rerunnable script with some added bells and whistles to check that files exist and ensuring any error causes the script to abort. These types of Bash scripts are quite simple to write: you’ve already learned important shell features like pipes, redirects, and background processes that play an important role in Bash scripts. In this section, we’ll cover the basics of writing and executing Bash scripts, paying particular attention to how create robust Bash scripts.

##### A robust Bash header

By convention, Bash scripts have the extension .sh. You can create them in your favorite text editor. Anytime you write a Bash script, you should use the following Bash script header, which sets some Bash options that lead to more robust scripts:

```
#!/bin/bash
set -e
set -u
set -o pipefail
```

These three options are the first layer of protection against Bash scripts with silent errors and unsafe behavior. Unfortunately, Bash is a fragile language, and we need to mind a few other oddities to use it safely in bioinformatics. We’ll see these as we learn more about the language.

##### Running Bash scripts

While we can run any script (as long as it has read permissions) with bash script.sh, calling the script as an executable requires that it have executable permissions. We can set these using:

```
$ chmod u+x script.sh
```

This adds executable permissions (+x) for the user who owns the file (u). Then, the script can be run with ./script.sh.

#### Variables and Command Arguments

Bash variables play an extremely important role in robust, reproducible Bash scripts. Processing pipelines having numerous settings that should be stored in variables (e.g., which directories to store results in, parameter values for commands, input files, etc.). Storing these settings in a variable defined at the top of the file makes adjusting settings and rerunning your pipelines much easier.

Unlike other programming languages, Bash’s variables don’t have data types. It’s helpful to think of Bash’s variables as strings (but that may behave differently depending on context). We can create a variable and assign it a value with (note that spaces matter when setting Bash variables—do not use spaces around the equals sign!):

```
results_dir="results/"
```

To access a variable’s value, we use a dollar sign in front of the variable’s name (e.g., $results_dir). You can experiment with this in a Bash script, or directly on the command line:

```
$ results_dir="results/"
$ echo $results_dir
results/
```

(399 to 423)

