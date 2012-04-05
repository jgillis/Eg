#!/bin/bash

buildbot create-master -r master
cp master/master.cfg.sample master/master.cfg
buildslave create-slave slave localhost:9989 example-slave pass

cd master && buildbot start && cd ..
cd slave && buildslave start && cd ..

