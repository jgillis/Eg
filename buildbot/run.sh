#!/bin/bash

rm -rf dummy dummy_rep master slave

echo "Creating a dummy project"
mkdir -p dummy.git && cd dummy.git && git --bare init   && cd ..
git clone dummy.git dummy

buildbot create-master -r master
cp master.cfg master/master.cfg

echo "tweak the settings a bit"

buildslave create-slave slave localhost:9989 example-slave pass

cd master && buildbot checkconfig master.cfg
cd ..
echo `pwd`
cd master && buildbot start && cd ..
cd slave && buildslave start && cd ..



firefox http://localhost:8010/

echo "You can go and login with dummy dummy"

cd dummy && echo "buildbot is awesome" > test.txt && git add test.txt && git commit -m "first commit" && git push origin master
echo "so far" >> test.txt && git add test.txt && git commit -m "second commit" && git push origin master && cd ..
