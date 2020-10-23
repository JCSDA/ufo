#!/bin/bash 

# clones specific branch of a given repo
# if branch does not exist clones develop
 
git_user=$1
git_token=$2
repo_name=$3
branch_name=$4
save_name=$5
save_dir=$6
branch_name_default=$7

echo "==============================================================================="
echo "Clone " $repo_name
echo "==============================================================================="

git ls-remote --heads --exit-code https://$git_user:$git_token@github.com/$repo_name $branch_name
exit_code=$?

if test "${exit_code}" == "0"; then
  branch_name_clone=${branch_name}
  echo ${branch_name} " branch found"
else
  branch_name_clone=${branch_name_default}
  echo ${branch_name} " branch does not exist"
  echo "clone " ${branch_name_clone}
fi

git clone -b $branch_name_clone https://$git_user:$git_token@github.com/$repo_name $save_dir/$save_name


