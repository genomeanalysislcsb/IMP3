## Bug report

<!-- Describe the bug/error you have encountered. -->

(Bug report)

### Log files and screenshots

<!-- Attach any relevant log files and screenshots, and their description. -->

(List of logs/screenshots)

## Steps to reproduce

<!-- How to reproduce the issue. -->

**Code version**
```txt
commit hash
branch name
```

**Config files**

<!-- Attach used and created config files. -->

(List of config files)

**Command**

<!-- Attach used launcher script and/or provide the command below. -->

```bash
# command used to launch IMP3
snakemake ...
```

**Input**

<!-- Provide relevant information about your input files. -->

(Used input)

**System**

<!-- Short description of the system setup where you are running `IMP3`. -->

(Used system)

## Possible fixes

<!-- If you can, link to the file or line of code that might be responsible for the problem - please make sure that the linked file corresponds to the indicated version of the code. -->

(Optional: How to fix)

/label ~bug

--------------------------------------------------

*NOTE: Please delete this section before submitting the issue.*

## Help

GitLab Markdown: [documentation](https://docs.gitlab.com/ee/user/markdown.html)

Formatting: please use [code spands and block](https://docs.gitlab.com/ee/user/markdown.html#code-spans-and-blocks) for
- shell commands and code
- console output, e.g. a log file snippet
- snippets from files, e.g. TXT, YAML, FASTA etc.

How to get the `IMP3` version:
```bash
# commit hash
git log -n 1 # the hash will be in the first line "commit ..."
# alternative option: human readable name
git describe --always

# branch name
git branch # the one with "*" is the current branch
```

Config files:
- any config files provided by the user to run `IMP3`
- file `sample.config.yaml` created by `IMP3` in the output folder
