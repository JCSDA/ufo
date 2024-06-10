#!/usr/bin/env bash
#
# (C) Copyright 2024 MetOffice Crown Copyright
#
# This script attempts to checkout a branch for a particular repository
# annotated in the body of the PR from which the GitHub action is triggered.
#
set -euo pipefail

action_repo="${GITHUB_REPOSITORY}"  # "owner/repository"
action_branch="${GITHUB_HEAD_REF}"
if [[ "${action_branch}" == 'develop' ]]; then
    exit
fi

target_dir="$1"
target_repo="${target_dir##*/}"
found=$(gh pr view "${action_branch}" -R "${action_repo}" \
        | grep 'build-group=' \
        | grep -E "/${target_repo}/|/${target_repo}#") || {
            echo "-- No build-group"; exit 0; }
echo "$found"

cd "$target_dir"
target_pr=${found##*/}
if [[ $target_pr =~ '#' ]]; then
    target_pr=${target_pr##*#}
fi
target_pr=${target_pr//$'\r'}  # remove those pesky invisible EOL characters
target_branch=$(gh pr list --search "$target_pr" \
                --json number,headRefName \
                -q ".[]|select(.number|contains($target_pr))|.headRefName")
git fetch origin "pull/${target_pr}/head:${target_branch}"
gh pr checkout "$target_pr"
echo "-- Switched to $(git branch --show-current)"

exit
