#!/usr/bin/env python3

import os
import sys
import stat
import tarfile
import urllib.request

bucket_name = "jedi-test-files"

download_file_name = sys.argv[1]
testfiles_name = sys.argv[2]
testfiles_path = sys.argv[3]
download_base_url = sys.argv[4]
md5check = sys.argv[5]

def DownloadUntar(download_base_url, download_file_name, testfiles_path, testfiles_name):
  urllib.request.urlretrieve( download_base_url+"/"+download_file_name+".md5", testfiles_path+"/"+testfiles_name+".md5")
  urllib.request.urlretrieve( download_base_url+"/"+download_file_name, testfiles_path+"/"+testfiles_name)
  tar_file = tarfile.open(testfiles_path+"/"+testfiles_name)
  tar_file.extractall(testfiles_path)
  tar_file.close()


if md5check == "1" :
  #  if .tar.gz and .tar.gz.md5 exist
  #  then download s3 md5
  #  and compare with local md5
  print ("checking md5sum")
  if os.path.isfile(testfiles_path+"/"+testfiles_name) and os.path.isfile(testfiles_path+"/"+testfiles_name+".md5") :
    print("local files found")

    #  dl md5 save it as *.md5.dl
    urllib.request.urlretrieve( download_base_url+"/"+download_file_name+".md5", testfiles_path+"/"+testfiles_name+".md5.dl")

    #  compare *md5.dl with md5 local
    with open(testfiles_path+"/"+testfiles_name+".md5", 'r') as f:
      md5_local = f.read()
    with open(testfiles_path+"/"+testfiles_name+".md5.dl", 'r') as f:
      md5_dl = f.read()
    if md5_local == md5_dl :
      print("no update in dataset")
    else:
      print("update found; download new dataset")
      DownloadUntar(download_base_url, download_file_name, testfiles_path, testfiles_name)
  else:
    print("local file not found; download from S3")
    print("downloading "+ download_base_url+"/"+download_file_name)
    DownloadUntar(download_base_url, download_file_name, testfiles_path, testfiles_name)

else:
  # downloading release data from DASH
  if os.path.isfile(testfiles_path+"/"+testfiles_name):
    print("local RELEASE file found")
  else:
    print ("dowloading RELEASE data from "+download_base_url+"/"+download_file_name)
    DownloadUntar(download_base_url, download_file_name, testfiles_path, testfiles_name)
