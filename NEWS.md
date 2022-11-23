# copyseparator 1.1.0

- Added a new parameter "rare_read" for copy_separate. During clustering analyses, clusters with less than this number of reads will be ignored. Default 10.
- If copy_separate found more than copy_number of clusters, this version will keep all major clusters (those have reads > rare_read). In the last version, copy_separate picks the copy_number of largest clusters.
- Sometimes some subsets may have no reads due to low coverage in the region. Those subset will now be skipped and the numbering of subsets will not be affected. Changes have also been made to copy_assemble to accommodate this change.
- When copy_number is 2 or 3, copy_assemble will identify the 2 or 3 most different clusters for those subsets that copy_separate retained > copy_number of clusters. When copy number >3, copy_assemble will still pick the copy_number of largest clusters.
- The input file for copy_assemble will be re-ordered to make sure the sequences are in the correct order (just in case that someone may have edited it).
- Remove the end part of a sequence if it contains too many ambiguous sites, which will otherwise affects the performance of copy_assemble.
- The shared gaps in the assembled final gene copy sequences will be removed.
- Added the new function sep_assem, which is a combination of copy_separate and copy_assemble. Now you can get the assembled final gene copy sequences by only running sep_assem. Only when the results are not satisfying, you can come back to make necessary changes and run coy_separate again.
- The input file and all resulting files will now be put in a single result folder. When you run multiple files one after another in a loop, the results will not be messed up.

# copyseparator 1.0.0
