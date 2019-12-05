ZZXAnalysis
==========

To install a complete CMSSW 10X area (including this package)
------------------------------
2016, 2017, and 2018 data analysis

Please use **CMSSW_10_2_18**. 

Download and execute the setup script with the following indications:
```
cmsrel CMSSW_10_2_18
wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/HHto4lX/ZZXAnalysis/master/checkout_10X.csh
cd $CMSSW_BASE/src
cmsenv
chmod u+x ${TMPDIR}/checkout_10X.csh
${TMPDIR}/checkout_10X.csh
scram b
```



To update this package from the release
------------------------------------------
In the package directory, simply issue
```
git pull
```

To commit and push new changes
------------------------------
To commit directly (you need write access to the repository):
```
git pull
[edit files]
```
Once you are ready to commit
```
git pull
git add [files to be added]
git commit -m ["commit message"]
git push origin master
```

Otherwise you can make a fork of the repository, develop therein, and make a pull request in the same way as for CMSSW.


