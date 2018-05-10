# **Welcome to the CDB lab**
## About settling in the linux environment for computational analyses

# YoSon Park  
June 6, 2015  

# overview

This post contains a quick overview of useful command line tools for research found in standard bash shells as well as some immediately useful tools. I'll try to go over some of very basic things for bash env in general and for linux cluster using LSF scheduler (which pmacs is). In general, example commands will be written out with comments (indicated with //).    
    
    

# pmacs setup quick ref

Some details can be found on the hpc website <http://www.med.upenn.edu/hpc/>  
Basic references for using pmacs including some LSF commands <https://hpcwiki.genomics.upenn.edu/index.php/HPC:User_Guide>  
I think pmacs lsf is v9. Full reference for lsf v9 is here <https://www-01.ibm.com/support/knowledgecenter/SSETD4_9.1.2/lsf_kc_cmd_ref.dita>  

# accessing pmacs from your local machine

From your local computer, you can use secure shell (ssh) to access pmacs. ssh is a cryptographic network protocol for initiating text-based shell sessions on remote machines. In principle, any network services can be securely accessed with ssh. The most relevant application in this lab is for accessing shell accounts on unix-like operating systems in a remote loation, e.g. pmacs.  


You will get to pmacs using this command  
ssh \<yourusername\>@consign.pmacs.upenn.edu  

## main login nodes (sometimes I'll call them headnodes)

File trasnfer processor node is separated within pmacs.  
mercury.pmacs.upenn.edu for all file transfer activities.  
consign.pmacs.upenn.edu for everything else.  


## notes about nodes

**Do not** run jobs on mercury.pmacs.upenn.edu unless it's a simple file or dir transfer  
**Do not** run jobs on consign.pmacs.upenn.edu unless you're accessing some LSF batch files unavailable on subnodes  



# tips on file or dir transfers to and from pmacs on your local machine
    
## rsync

If you need to move over any files from your computer to pmacs, I recommend you use rsync. There are other methods with pros and cons. rsync is useful because it combines file synchronization as well as transfer within unix-like systems. Also it generally won't leave traces of partially transferred/synced files. So if the connection was interrupted for any reason, it will terminate and remove the partially transferred file.  
    
The basic usage is that you prompt rsync, then indicate origin and destination. paths ending with / will assume the destination is final. paths not ending with / will create dir mirroring the last branch of the dir tree to the destination.  

## rsync example

//copies my\_important\_file.txt to your home dir on pmacs. -a indicates that it will preserve all metainfo of files, e.g. permission levels and timestamps  
rsync -a my\_important\_file.txt \<yourusername\>@mercury.pmacs.upenn.edu:/home/\<yourusername\>/  
    
//mirrors the local workdir directory with your home dir on pmacs by creating workdir if non-existent and synchronizes files (overwrite if source newer) if dir already existant.  
rsync -a workdir \<yourusername\>@mercury.pmacs.upenn.edu:/home/\<yourusername\>  
    
You could do the reverse by simply reversing the order of two arguments from above. Alternatively, if you're already on pmacs and want to do some file transfers to your local machine, you need to know your ip. ip is a code that functions as an address of your system, e.g. computer, in the network. There are different versions but the usual one you'd see is 4 sets of 1-3 digits of numbers separated by a period (.) called IPv4, e.g. 130.91.#.###  

pmacs does not have a file backup system. 
**Do not** overwrite files others may be using. Be careful and make sure when you are synchronizing dirs or files.  


# nodes other than consign / mercury

Each "node" is like a computer. It has memory cards and CPUs (sometimes GPUs). There are suposedly ~5k virtual cores in pmacs. These cores are organized into ~150 nodes (consign and mercury being two of them and most others are named with numeric ids) and each contain 24-64 cores. This organization is irrelevant to you unless you're running jobs for multiple-threads. These nodes are then organized into queues, mainly normal and denovo. You probably will exclusively be using normal queue. You'll hear me rambling about cores and nodes and threads interchangeably sometimes. It suffices to say you won't have to worry about waiting for your job to be submitted for too long. Instead of having its own hard drive like your computer, it's connected to a file system, allowing your files to be accessible no matter which node you end up using (well, most of the time anyway).  


## accessing interactive nodes

So the first thing once you're setup with a pmacs account is probably connecting to the subnode where you can do stuff.

bsub -Ip bash  
bsub -Is bash  

Either of the above commands will "submit" a job to consign, which in this case simply means you've relocated to the node from consign. Now you are on one of what's called interactive nodes of pmacs. I think there's a limit set at 10 interactive nodes per user. The rule of thumb is that if you have to run things with 10 open terminals at the same time, you probably should write a script to submit it to other nodes.  

## accessing other nodes

For all nodes that will do heavy lifting for you, you may write a script and send it over using bsub.  
example script bogus.sh  
#!/bin/bash  
echo "hello cdb lab!"  

then you may submit it as such    
bsub -e bogus.err -o bogus.out -M 10 -W 1:00 < bogus.sh  

for bsub,   
-e // $STDERR redirect file  
-o // $STDOUT redirect file  
-M // memory limit in Mb  
-W // hour limit in the form of HH:MM  


For other parameters, check on pmacs  

man bsub  

Any job without -q parameter will be submitted to the normal queue by default on pmacs. Default -M for norml queue is set at 6000 Mb (6Gb). Depending on the job, you may adjust it up to 250Gb. For instance, star for rna-seq alignment will use 30 - 50Gb memory.  


# software installation and usage in general 

Linux clusters use what's called modules to organize different versions and types of software for users. Often a well-organized cluster will only have super archaic and can't-live-without-it type tools installed globally. For everything else, either you install it yourself in your home dir or you find a module that loads a virtual environment connecting you to it.  

You can find modularly installed software on pmacs here <https://hpcwiki.genomics.upenn.edu/index.php/HPC:Software>  

For tools that not many are using, you can attempt to contact pmacs hpc admin to install it for you. Likely, you'll have to install it yourself. Depending on how frequently we all would be using it, I may install it for everyone on the project space as necessary. Otherwise, use /home/\<yourusername\>/ and locally install what you need for your projects.  
    

## looking for docs  

Generally, each software dir, for example, should contain a doc such as README or README.md that explains some basics about what's in the dir. README is a text file you may open using vi, cat, less, more, or your favorite commandline text editor. README.md is also a text file but sometimes it's a doc like the one you're reading right now that has been generated using a markdown editor. For most purposes, it can be accessed just like any text file. You may see some weird (non-english) lines dictating typesettings, etc.  
    
More globally, you may use commands such as man and info to see whether the tool has a text version of the manual. This is handy particularly when you are using some older tools that you may not know specific versions of. Google search works for general subjects but what's on the system man or info would be the most relevant doc for what you're using as often parameters can be deprecated, etc.  
    
example:  
man bsub  
info cut  
    

# working with modules

While on pmacs, you may list all available installed modules as such  
module avail

For example, if you don't see R 3.1.2 anywhere and you really want to use R 3.1.2 and no other versions, then you'd approach it by  
// check existing global installation of R  
R --version  
// see if your favorite version is available on pmacs as a module  
module avail  
// since it's available, you may check any additional information regarding this module installation  
module show R-3.1.2  
// now that you know that there's no specific strings attached to this module, you may load it to start using it  
module load R-3.1.2  
// check to see if you've successfully loaded R-3.1.2  
R --version  

The same can be done for python, or for specific tools very useful to us such as star and sailfish, or their related libraries such as zlib, boost, etc.  

# automatically loading modules

If you're like me and can't spend a day without checking something using python and you get annoyed every time you realize you forgot to load it, you may add this line to your .bashrc file. .bashrc if a configuration file used by the operating system and often used with .profile and .bash_profile files to store your preferences (or lack thereof). On pmacs, you may see if you have one in /home/\<yourusername\>/ and you can create one if you don't. Configuration files are generally loaded in a location or context specific manner.  

if [ $HOSTNAME != "consign.hpc.local" ] && [ $HOSTNAME != "mercury.pmacs.upenn.edu" ]; then  
  module load python-2.7.9  
fi  

Then the python-2.7.9 module will be automatically loaded whenver you're on pmacs  

    
# miscellaneous info about working on linux

In linux environment, angle brackets direct how data flows.  
">" indicates what's on the left is passed on to what's on the right.  
// cat displays the content of input.txt (referred to as $STDOUT in linux env) which then is stored as same.input.txt  
cat input.txt > same.input.txt  
"<" for the opposite effect.  
//submit a job using myscript.sh to the pmacs cluster  
bsub < myscript.sh  


# rules on working in the shared space  

This may sound very basic and obvious to you. I cannot stress enough how important it is to make sure shared dir remain safe and organized for everyone's sanity (especially mine).  
    
Much of mistakes attributable to corrupted input files, wrong versions of input files, etc. can be prevented or detected early simply by having a good naming convention for files and dirs. Files we use daily can range from your favorite shell script to genomic annotations to sequencing data files. For some of these files, using a different version, for instance, is not very obvious to figure out. Generally speaking, if you don't know what the file is and you don't know whether anybody's using it, **do not** relocate or overwrite it. If you think you might forget what it is for or someone else may be using it, leave some notes in the dir to indicate what it is. It's o.k. to make mistakes like this but make it convinient to find out when and where you did. Try to avoid common names like dbsnp.txt when you really meant dbsnp.build142.txt  


# awk & sed

awk (also implemented in nawk and gawk) and sed allow easy manipulations of text files. awk is often used for delimited files. sed works better with streams of texts. The syntaxes are relatively crude and much simpler than perl or python but still very fast. awk is extremely useful for operating on columns of data.  
intermediate awk tutorial can be found here <http://www.grymoire.com/Unix/Awk.html>  

## toy examples

// print the second and third columns separated by a tab from input.txt and redirect the results in output.txt

awk '{print $2, "\\t", $3}' input.txt > output.txt

// find values that are negative and replace each column in input.txt with 0 and store results in output.txt

awk '{for (i=1; i<=NF; i++) if ($i < 0) $i = 0; print}' input.txt > output.txt

// list the number of columns in each row of input.txt

awk '{print NF}' input.txt

// substitutions

awk '{print toupper(\$0)}' input.txt > output.txt // changes lower cases to upper cases  
awk '{print tolower(\$0)}' input.txt > output.txt // changes upper case to lower cases  
sed 's/cat/dog/g' input.txt > output.txt //substituting "cat" in input.txt to "dog"  
awk '/animal/ { gsub(/cat/, "dog") }; { print }' input.txt > output.txt // same as above only for lines containing word "animal"  
sed 's/[ \\t]\*$//' input.txt > output.txt // delete all trailing white space and tab in a file.  
sed '/^\\$/d' input.txt // delete blank lines in input.txt  

# working with patterns and numbers in a sequence file

// sometimes pipe (|) is used to bypass having to write out an intermediate file in linux

zcat input.vcf.gz | awk '($1!\~"##"){print $1, $2}' > output.txt //skip vcf headers and print out chromosome and base position from a compressed vcf file  
awk '($1=="chr1" && $2>1000 && $2<5000){print}' input.vcf // print all rows in vcf file if the variant position is between 1000 and 5000 on chromosome 1  
awk '/ATGCATGCATGC/ { n++ }; END { print n+0 }' input.txt > output.txt // pick out the total number of lines containing the pattern ATGCATGCATGC  
sed -n '/^.\{76\}/p' input.txt // print lines that are longer than 76 characters  
sed -n '1001p' input.txt // print only line number 1001  
sed -n '950,1049p' input.txt // print lines 950 to 1049 (inclusive)  
awk '{sum+=$9} END {print sum}' input.txt // sum column 9 of the input.txt  
// print total number of reads, total number unique reads, percentage of unique reads, most abundant sequence, its frequency and percentage of total in input.fastq.gz  
zcat input.fastq.gz | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique\*100/total,maxRead,count[maxRead],count[maxRead]\*100/total}'  
sed -n '2\~4p' input.fastq  // print every 4th line starting at the second line (extract the sequence from input.fastq)  


# other useful oneliners

wc -l // count number of lines in a file  
ls | wc -l // count number of files in a directory  
tac input.txt // print the file in reverse order  
sed 's/^M$//' input.txt // converts a dos file into unix mode  
nl input.txt > numbered.input.txt // add a column counting line numbes of the file
sed = input.txt | sed 'N;s/\\n/ /' > numberd.input.txt // same as above in a longer command  
shuf input.txt | head -n 10 // print 10 random lines from input.txt  
echo {A,T,G,C}{A,T,G,C}{A,T,G,C} // print all possible 3mer dna sequence combinations  
find . -name "\*.bam" // search for \*.bam files in the current dir recursively  



