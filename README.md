# tfProfiling

tfProfiling is a group of convenient tools based on [TranscriptionFactorProfiling](https://github.com/PeterUlz/TranscriptionFactorProfiling) published by PeterUlz for the paper [Inference of transcription factor binding from cell-free DNA enables tumor subtype prediction and early detection](https://www.nature.com/articles/s41467-019-12714-4) that further improves the experience of analyzing Transcription Factors (TFs) in cell-free DNA (cfDNA) using R.

## 1. Description

tfProfiling enhances the process of analyzing TFs profiles in cfDNA data by providing a simplified installation process for the required bioinformatics tools, as well as providing a series of functions that allow calling these tools in R coding language. tfProfiling can be used, either as a standalone tool or in combination of other tools, such as [ULPwgs](https://github.com/TearsWillFall/ULPwgs) or [DNAfrags](https://github.com/TearsWillFall/DNAfrags) to create complex workflows to analyze cfDNA data.



![Workflow for tfProfiling](https://github.com/TearsWillFall/tfProfiling/blob/master/Graph.png?raw=true)


**Tools:**
* [bedtools](https://github.com/arq5x/bedtools2): Bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks
* [htslib](https://github.com/samtools/htslib): C library for high-throughput sequencing data formats 
* [samtools](https://github.com/samtools/samtools/):Tools (written in C using htslib) for manipulating next-generation sequencing data


## 2. System Requirements

Tested on fresh Ubuntu and Arch Linux installations. Should be working on other Linux distros too as long as equivalent packages are provided. 

In order to be able to download and compile the source files of all the required tools the following programs are required:

* git
* make
* gcc
* cmake
* autoconf

These tools can and should be installed using the terminal with the following commands:


* **For Ubuntu:**

  ```
  sudo apt install git make gcc ant cmake autoconf 
  ```

* **For Arch Linux:**

  ```
  sudo pacman -S git ant make cmake gcc autoconf
  ```
  
Additional dependencies may be needed to succesfully install `devtools` package in R:
  
* **For Ubuntu:**

  ```
  sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
  ```
  
* **For Arch Linux:**

  ```
  sudo pacman -S build-essential libcurl-gnutls libxml2-dev openssl
  ```
  
  
 ## 3. Installation Instructions

In order to install `tfProfiling` package we will be using R `devtools`:

```
install.packages("devtools")
devtools::install_github("TearsWillFall/tfProfiling")
```
If `devtools` package installation fails check System Requirements section, as you may be missing a dependency.

Once the `tfProfiling` package is installed we can use the function `install_required_tools()` to download and set up all the tools required for the bioinformatic process. This will create a directory named `tools` in the current working directory with all the tools. **Note: All functions within this package call scripts from the `tools` directory, therefore if this directory is moved, deleted or the current working directory is changed, these function will fail.**

```
tfProfiling::install_required_tools()
```

Or, alternatively.

```
library("tfProfiling")
install_required_tools()
```

## 4. Usage

![Workflow for tfProfiling](https://github.com/TearsWillFall/tfProfiling/blob/master/Graph2.png?raw=true)


