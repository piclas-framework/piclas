# GitLab Workflow

Code development is performed on the [GitLab platform](https://gitlab.com/piclas/piclas), with the protected `master` and
`master.dev` branches. The actual development is performed on feature branches, which can be merged to `master.dev` following a
merge request and the completion of a merge request checklist. After a successful pass of the nightly and weekly regression test,
the `master.dev` can be merged into the `master`. A merge of the `master.dev` to the `master` should be associated with a release
tag, where the changes to previous version are summarized.

In the following the envisioned development process using issues and milestones, the release & deploy procedure as well as other
developer relevant issues are discussed.

## Issues & Milestones
Issues are created for bugs, improvements, features, regression testing and documentation. The issues should be named with a few
keywords. Try to avoid describing the complete issue already in the title. The issue can be assigned to a certain milestone
(if appropriate).

Milestones are created based on planned releases (e.g. Release 1.2.1) or as a grouping of multiple related issues (e.g.
Documentation Version 1, Clean-up Emission Routines). A deadline can be given if applicable. The milestone should contain a short
summary of the work performed (bullet-points) as its contents will be added to the description of the releases. Generally, merge
requests should be associated with a milestone containing a release tag, while issues should be associated with the grouping milestones.

As soon as a developer wants to start working on an issue, she/he shall assign himself to the issue and a branch and merge request
denoted as work in progress (`WIP: ...`) should be created to allow others to contribute and track the progress. For this purpose,
it should be created directly from the web interface within the issue (`Create merge request`). This creates a branch, naming it
automatically, leading with the issue number (e.g. `60-fix-boundary-condition`) and associates the branch and merge request to the
issue (visible in the web interface below the description). To start working on the issue, the branch can be checked out as usually.

Ideally, issues should be created for every code development for documentation purposes. Branches without an issue should be avoided
to reduce the number of orphaned/stale branches. However, if branches are created outside of the issue context, they should be named
with a prefix indicating the purpose of the branch, according to the existing labels in GitLab. Examples are given below:

    feature.chemistry.polyatomic
    improvement.tracking.curved
    bug.compiler.warnings
    reggie.chemistry.reservoir
    documentation.pic.maxwell

Progress tracking, documentation and collaboration on the online platform can be enabled through creating a merge request with the
WIP prefix for this branch instead of an issue. An issues created afterwards cannot be associated with an already created branch,
without renaming the branch to include the issue number at the beginning. However, this should be avoided.

## Merge Request

Merge requests that are not WIP are discussed every Monday by the developer group to be considered for a merge. The following
checklist has to be completed before a merge request should be approved. For bugs only the first points have to be considered,
while for features and improvements the complete list has to be completed.

* [ ] Style Guide
* [ ] Maximum of 10 compile warnings via *./tools/test_max_warnings.sh*
* [ ] Descriptions for new/changed routines
  * Short header description (do not just spell out the name of the subroutine, units for important variables if applicable)
  * Workflow (short summary in the header, inside the routine at the appropriate positions)
* [ ] Reggie (small test setup, entry in REGGIE.md table, readme.md in test case folder)
* [ ] New feature description in appropriate documentation (user/developer guide)

For this purpose, the developer can select the respective template for his merge request (Bug: only first two to-do's or Feature:
all five to-do's, Improvements can utilize either depending on the nature of the improvement). The appropriate checklist will then
be displayed as the merge request description. Merge requests generated automatically through the Issues interface have already
`Closes #55` as a description. When editing the merge request, the description gets overwritten by the template. Thus, the issue
number has to be added manually after the template is chosen. The templates for merge requests are stored in *./.gitlab/merge_request_templates/*.

## Release and deploy

A new release version of PICLas is created from the **master** repository, which requires a merge of the current **master.dev**
branch into the **master** branch. The corresponding merge request should be associated with a release milestone (e.g. *Release
1.X.X* within which the merge request is referenced, e.g., "See merge request !283"). Within this specific merge request, the
template `Release` is chosen, which contains the to-do list as well as the template for the release notes as given below.
After the successful completion of all to-do's and regression checks (check-in, nightly, weekly), the **master.dev** branch can be
merged into the **master** branch.

### Release Tag

A new release tag can be created through the web interface ([Repository -> Tags](https://gitlab.com/piclas/piclas/tags) -> New tag)
and as the `Tag name`, the new version number is used, e.g., 

    v1.X.X

The tag is then created from the **master** branch repository and the `Message` is left empty. The release notes, which were used 
within the corresponding milestone, shall be given in the following format

    ## Release 1.X.X

    ### Documentation

    * Added section about particle emission

    ### Reggie

    * Added a regression test of the chemistry routine

    ### Features

    * Skipping field update for the HDG solver for user-defined number of iterations

    ### Improvements

    * Speed-up by skipping/cycle over neutral particles in deposition

    ### Fixes

    * Treatment of non-linear polyatomic molecules during analyze and wall interaction

Headlines without changes/additions within a release can be omitted.

### Collaborative Numerics Group

The collaborative numerics group for PICLas is located at *https://gitlab.com/collaborative-numerics-group/piclas*.
There are two possible ways of sharing code with other people and are explained in the following.

#### Single master branch repository
Method 1 involves a single master branch that is updated from the internal PICLas repository, similar
to the repository of PICLas at *github.com*.

The master branch of development group can be merged after the successful regression check with the master of the collaborative group.
For this purpose, the collaborative repository can be added as a remote (this step has only to be performed once)

    git remote add remote_name git@gitlab.com:collaborative-numerics-group/piclas/piclas.git

First, make sure to have the most recent version of the master branch (of the development repository)

    git checkout master && git pull

Now you can checkout the most recent version of the master branch of the collaborative-numerics-group and create a local branch
with that version (performing only a simple checkout will create a detached HEAD state)

    git fetch
    git checkout -b branch_name remote_name/master

The master branch of the development repository can now be merged into the newly created branch

    git merge origin/master

Finally, the changes can be pushed from the local branch *branch_name* to the master of collaborative-numerics-group

    git push remote_name master

If a tag has also been created, it should be pushed separately.

    git push remote_name tag_name

Afterwards, the local branch *branch_name* can either be deleted or utilized for future merges

    git branch -d branch_name

#### Pushing each branch to a separate repository
Method 2 involves pushing each development branch from the internal PICLas repository to a separate
repository in the PICLas CRG separately.

Extract solely a single branch by cloning only *myfeaturebranch* via

    git clone -b myfeaturebranch --single-branch git@gitlab.com:piclas/piclas.git piclas_new_CRG ;

Navigate to the newly created directory

    cd piclas_new_CRG

Push this repository that only contains one branch to the CRG via

    git push --mirror git@gitlab.com:collaborative-numerics-group/piclas/piclas.myfeaturebranch.git

which created a new repository that only contains the desired feature branch. When sharing this new
repository with new people, simply add them as members of this new repository (not the complete
CRG!).


### GitHub

Finally, the release tag can be deployed to GitHub. This can be achieved by running the `Deploy` script in the
[CI/CD -> Schedules](https://gitlab.com/piclas/piclas/pipeline_schedules) web interface. At the moment, the respective tag and the
release have to be created manually on GitHub through the web interface with **piclas-framework** account. The releases are
accessed through [Releases](https://github.com/piclas-framework/piclas/releases) and a new release (including the tag) can be
created with `Draft a new release`. The tag version should be set as before (`v1.X.X`) and the release title accordingly
(`Release 1.X.X`). The release notes can be copied from the GitLab release while omitting the `## Release 1.X.X` headline as it was
given with the release title before.
