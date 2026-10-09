# Development Workflow

Code development is performed on the [GitLab platform](https://piclas.boltzplatz.eu/piclas/piclas), with the protected `master` and
`master.dev` branches. The actual development is performed on feature branches, which can be merged to `master.dev` following a
merge request and the completion of a merge request checklist. After a successful pass of the nightly and weekly regression test,
the `master.dev` can be merged into the `master`. A merge of the `master.dev` to the `master` should be associated with a release
tag, where the changes to previous version are summarized.

In the following, the envisioned development process using issues, milestones and merge requests is discussed.
The release procedure is described in Chapter {ref}`developerguide/release:Release & Deployment`.

## Issues & Milestones
Issues are created for bugs, improvements, features, regression testing and documentation. The issues should be named with a few
keywords. Try to avoid describing the complete issue already in the title. The issue can be assigned to a certain milestone
(if appropriate).

Milestones are created based on planned releases (e.g. Release 1.2.1) or as a grouping of multiple related issues (e.g.
Documentation Version 1, Clean-up Emission Routines). A deadline can be given if applicable. Generally, at least merge
requests should be associated with a release milestone.

As soon as a developer wants to start working on an issue, she/he shall assign himself to the issue and a branch and merge request
denoted as work in progress (`Draft: ...`) should be created to allow others to contribute and track the progress. For this purpose,
it should be created directly from the web interface within the issue (`Create merge request`). This creates a branch, naming it
automatically, leading with the issue number (e.g. `60-fix-boundary-condition`) and associates the branch and merge request to the
issue (visible in the web interface below the description). To start working on the issue, the branch can be checked out as usually.

If branches are created outside of the issue context, they should be named with a prefix indicating the purpose of the branch,
according to the existing labels in GitLab. Examples are given below:

    feature.chemistry.polyatomic
    improvement.tracking.curved
    bug.compiler.warnings
    reggie.chemistry.reservoir
    documentation.pic.maxwell

Progress tracking, documentation and collaboration on the online platform can be enabled through creating a merge request with the
`Draft:` prefix for this branch instead of an issue.

## Merge Request (MR)

A merge request can be created to track the development of the branch with the intent to merge the development.
Merge requests that are not WIP are discussed every Monday by the developer group to be considered for a merge.
In prepartion of the merge, the developer can select the respective template for the merge request (**Bug** or **Feature**, stored in `.gitlab/merge_request_templates/`).
Improvements can utilize either depending on the nature of the improvement. The appropriate checklist will then
be displayed as the merge request description. Note that merge requests generated automatically through the Issues interface have already
`Closes #55` as a description. When selecting a template, the description gets overwritten and the issue
number has to be added back manually.

The merge request title should be of the form

    [feature.branch.name] Short and concise feature/improvement/bugfix description

and will be automatically utilized within the release notes based on the selected label for the merge request. MR labelled `Feature`, `Improvement`, and `Bug` will be automatically sorted under the respective release notes headings. Every other label will be included in the `Documentation/Tools/Regression Testing` section. To avoid the addition of a merge request to the release notes (e.g. repetitive merge requests testing the Gitlab CI or similar), use the `No-release-notes` label.

### Automatic release note generation

The Gitlab CI configuration `.gitlab-ci.yml` includes the `update_release_notes` job, which runs for every feature branch merge into `master.dev`. The job determines the merge request belonging to the current pipeline (which provides the commit SHA) and extracts the merge request title and label to add to a merge request with the `Release` label. If such a merge request is not available, it will be created. The contents of the job can be tested locally by setting the following variables

    export CI_API_V4_URL="https://piclas.boltzplatz.eu/api/v4"
    export CI_PROJECT_ID="X"                # Can be provided upon request
    export RELEASE_BOT_TOKEN="glpat-XXX"    # Can be provided upon request
    export CI_COMMIT_SHA="4e39ce852"        # Select a recent merge commit

and copying everything within `script:` section below `set -e` (not including it) to a separate script file (e.g. `update_release_notes.sh`). Make sure to comment out the last `curl -sf --request PUT` command, when testing to avoid modifying the description of a current merge request:

    # Update the release MR description
    # curl -sf --request PUT \
    # --header "PRIVATE-TOKEN: ${RELEASE_BOT_TOKEN}" \
    # --header "Content-Type: application/json" \
    # --data "$(jq -n --arg d "$UPDATED" '{description: $d}')" \
    # "${API}/merge_requests/${REL_IID}" > /dev/null

    # Instead of updating, just show what would happen
    echo "DRY RUN - would add to !${REL_IID}:"
    echo "$ENTRY"
    echo ""
    echo "Full updated description:"
    echo "$UPDATED"

The script assumes that the default Gitlab merge request commit message containing `See merge request !123` has not been modified.
