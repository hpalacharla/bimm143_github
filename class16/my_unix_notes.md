---
output: html_document
editor_options: 
  chunk_output_type: inline
---

Nano is just a basic text editor

Core unix commands include:

- pwd: Print working directory
- ls: List files in the current directory
- cd: Change directory
- mkdir: Make a new directory (folder)
- touch: Create a new file
- rm: Remove/Delete a file
- nano: A basic text file
- mv: Move/Rename a file
- cp: Copy/Duplicate a file
- less: View a file (Q to quit)
- man
- head/tail

Key point: pwd returns the "absolute path", where you are relative to the root of the file system, eg: `/Users/harshi/Desktop/class16`

Useful shortcuts include

- ../ Parent dir
- ~/  Home dir
- *   Wild card match everything 

## AWS instance

ssh -i ~/Downloads/bimm143_harshitha.pem ubuntu@ec2-16-146-121-254.us-west-2.compute.amazonaws.com

scp -i ~/Downloads/bimm143_harshitha.pem ubuntu@ec2-16-146-121-254.us-west-2.compute.amazonaws.com:~/work/myresults.txt .

scp -i ~/Downloads/bimm143_harshitha.pem ubuntu@ec2-16-146-121-254.us-west-2.compute.amazonaws.com:~/work/mm-third.x.zebrafish.txt .

scp -i ~/Downloads/bimm143_harshitha.pem ubuntu@ec2-16-146-121-254.us-west-2.compute.amazonaws.com:~/work/mm-third.x.zebrafish.tsv .

# Make sure I am in my class16 dir
