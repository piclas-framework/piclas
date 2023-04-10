## Related Issue

Closes #number

## Merge Request Checklist

* [ ] Style Guide
* [ ] Maximum of 10 compile warnings via *./tools/test_max_warnings.sh*. How many warning were found?
* [ ] No large files via *./tools/test_max_file_size.sh*. What is the largest file?
* [ ] Descriptions for new/changed routines
  * [ ] Short header description (do not just spell out the name of the subroutine, units for important variables if applicable)
  * [ ] Workflow
    * [ ] Short summary in the header
    * [ ] Inside the routine at the appropriate positions
* [ ] Reggie
  * [ ] Add small test setup
  * [ ] Add entry in REGGIE.md table
  * [ ] Check automatic restart functionality of reggie example via Load Balance (checks correct allocation and deallocation for the test case)
* [ ] New feature description in appropriate documentation (user/developer guide)
* [ ] Check that no large files were added to the repository
