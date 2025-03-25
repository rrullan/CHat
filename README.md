# CHat
Characterization and Analysis of electronic Transition software

The obective of these scripts is to compute some indexes to help us characterize the transitions both for solid and molecular systems.
For now, only the Tozer index has been thoroughly studied and tested. The DCT index CP2K implementation is still work to do. 
We will add the tools as time goes. 

How to use it: 
 -Having the excited\_index.sh and excited\_index.py in the directory.
 -For CP2K: you should have the file.inp, file.out and the mofile.cube in the directory.
 -For Orca: you should have the file.out, file.gbw and eventually the file.cis in the directory.
            It should work either for Orca 5 and Orca 6.

WARNING: for Orca, you should give the path to orca\_plot for your architecture in excited\_index.sh

Once you have all the files:
bash excited\_index.sh 

WARNING: you should have the extension in the name of the files 


This code uses the Cube-Toolz package to manipulate cube files (see the github repository: https://github.com/funkymunkycool/Cube-Toolz).
This code is under the license: European Union Public License 1.2.

